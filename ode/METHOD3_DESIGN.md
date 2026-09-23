# Method 3: Optical Depth Gradient Monitoring

## Design Specification

### Overview
Method 3 monitors the optical depth gradient dτ/dl during integration and dynamically constrains the step size when the gradient is large, preventing the adaptive stepper from taking steps that skip over resonance layers.

### Implementation Details

#### 1. Access to dτ/dl
The optical depth gradient dτ/dl is available as `dydx[3]` (index 3 in the state derivative vector), computed by `PhotonEvolution::calc_dtaudl()`. The Dormand-Prince 5(4) stepper evaluates the derivative function at 7 stages per step:
- `dydx_old` (start of step)
- `k2, k3, k4, k5, k6` (intermediate stages)
- `dydx_new` (end of step)

#### 2. Gradient Threshold and Limits
**Parameters:**
- `dtau_threshold`: Threshold value for "large" dτ/dl (default: 10.0 km⁻¹)
- `max_dtau_per_step`: Maximum allowed Δτ accumulation per step (default: 0.1)
- `tau_index`: Index of τ in state vector (default: 3)

**Logic:**
```
max_dtaudl = max(|dydx_old[3]|, |k2[3]|, |k3[3]|, |k4[3]|, |k5[3]|, |k6[3]|, |dydx_new[3]|)

if max_dtaudl > dtau_threshold:
    h_gradient_limit = max_dtau_per_step / max_dtaudl
    h_proposed = min(h_proposed_by_PI_controller, h_gradient_limit)
```

#### 3. Step Size Adjustment Formula
When dτ/dl exceeds the threshold, the step size is constrained such that:
```
h < max_dtau_per_step / max_dtaudl
```

This ensures:
```
Δτ ≈ h × dτ/dl < max_dtau_per_step
```

#### 4. Integration with PI Controller
The gradient-based limit is applied **after** the standard PI (Proportional-Integral) controller computes its proposed step size. The final step size is:
```
h_new = min(h_PI_controller, h_gradient_limit)
```

This ensures both error control (from PI controller) and gradient control (from resonance monitoring) are satisfied.

#### 5. Stage Checking
Method 3 checks dτ/dl at **all 7 Runge-Kutta stages**, not just endpoints. This is critical because:
- A resonance layer might be encountered mid-step
- Using only endpoint values could miss narrow peaks in dτ/dl
- The intermediate stages (k2-k6) sample the derivative at different points within the step

### Pseudocode

```cpp
bool success(double err) {
    // Standard PI controller
    if (err <= 1.0) {
        double scale = compute_PI_scale(err, errold);
        double h_proposed = h_old * scale;
        
        // Gradient monitoring
        double max_dtaudl = max({
            abs(dydx_old[tau_index]),
            abs(k2[tau_index]), abs(k3[tau_index]),
            abs(k4[tau_index]), abs(k5[tau_index]),
            abs(k6[tau_index]), abs(dydx_new[tau_index])
        });
        
        // Apply gradient constraint if needed
        if (max_dtaudl > dtau_threshold) {
            double h_gradient_limit = max_dtau_per_step / max_dtaudl;
            h_proposed = min(h_proposed, h_gradient_limit);
        }
        
        h_new = h_proposed;
        return true;  // Step accepted
    } else {
        // Step rejected, reduce h_old
        return false;
    }
}
```

### Theoretical Justification

Near a cyclotron resonance (ω = ω_c), the scattering cross-section and hence dτ/dl can become very large. If the adaptive stepper takes a large step h across the resonance, it will significantly underestimate the accumulated optical depth:

```
τ_true = ∫ (dτ/dl) dl ≈ ∫ peak_value dl
τ_computed ≈ h × average_value << τ_true
```

By constraining h when dτ/dl is large, we ensure adequate sampling of the resonance layer:
```
h < ε / (dτ/dl)_max
```
where ε = max_dtau_per_step is the tolerance for optical depth accumulation per step.

### Parameter Selection Guidelines

**dtau_threshold:**
- Too low: Unnecessarily constrains step size everywhere, reducing efficiency
- Too high: Fails to detect resonances, allowing large steps across them
- Recommended: 1.0 - 10.0 km⁻¹ (depends on typical background dτ/dl)

**max_dtau_per_step:**
- Controls the resolution of optical depth accumulation
- Smaller values: More accurate but more steps (slower)
- Larger values: Faster but may miss narrow resonances
- Recommended: 0.01 - 0.1 (since target scattering optical depths are O(1))

### Current Benchmark Results

Initial tests show that with the default photon initialization (init07), the maximum encountered dτ/dl is very small (< 1 km⁻¹), well below the threshold. This indicates:

1. **Either:** The photons in these test cases don't encounter strong resonances
2. **Or:** The resonances are being missed entirely (the problem we're trying to solve)

To properly validate Method 3, we need test cases where photons definitely pass through resonance layers with high dτ/dl.

### Comparison with Other Methods

**Method 1 (Original - Events only):**
- Relies only on event detection (when τ crosses threshold)
- Can take arbitrarily large steps between events
- May miss resonances if step crosses entire layer

**Method 2 (h_max constraint):**
- Simple: limits h ≤ h_max always
- Non-adaptive: same constraint everywhere
- Inefficient: many small steps even in smooth regions
- For h_max = 0.01 km, takes ~1M steps vs ~50 steps for Method 1

**Method 3 (Gradient monitoring):**
- Adaptive: only constrains h when dτ/dl is large
- Physics-aware: directly monitors the quantity that matters
- Efficient: allows large steps in smooth regions, small steps near resonances
- In current tests: Same efficiency as Method 1 (gradient never triggered)

### Implementation Files

- **Stepper:** `/Users/tz_mbp/Desktop/magmc_tz/ode/dopr5_gradient.hpp`
- **Benchmark:** `/Users/tz_mbp/Desktop/magmc_tz/ode/benchmark.cpp`
- **Detailed benchmark:** `/Users/tz_mbp/Desktop/magmc_tz/ode/benchmark_detailed.cpp`

### Next Steps

To properly test Method 3, we need:

1. Identify photon initial conditions that produce high dτ/dl (near resonances)
2. Adjust threshold parameters based on actual dτ/dl scales in the problem
3. Compare computed scattering optical depths across methods for accuracy
4. Verify that Method 3 catches resonances that Method 1 misses

### Usage Example

```cpp
// Create stepper with gradient monitoring
StepperDopr5Gradient<4, PhotonEvolution> stepper(
    photon_evolution,
    1e-6,        // atol
    1e-6,        // rtol
    0.1,         // max_dtau_per_step
    10.0,        // dtau_threshold
    3            // tau_index
);

// Use like normal stepper
stepper.init(x0, h0, y0);
while (...) {
    stepper.do_step();
    int event = stepper.detect_event();
    // ... handle events
    stepper.update_old();
}
```
