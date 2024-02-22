# LU-Factorization-with-BANs

The utilization of LU factorization aimed to address linear systems involving bounded arithmetic and both infinite and infinitesimal perturbations. However, the test results indicate that the algorithm introduces significant numerical errors in specific scenarios, and the Bounded Arithmetic Numbers (BANs) exhibit challenges in achieving optimal solutions.

🔍 **Issues Identified:**
- Problematic cases were identified in scenarios with an exceptionally high condition number and infinitesimal coefficients during factorization.
- Instances of a simplex loop occurred, causing the algorithm to fail in converging to any solution in its initial implementation.
- The introduction of a threshold helped alleviate issues in specific cases but proved less effective when dealing with infinite coefficients.

🛠️ **Algorithm Modifications:**
- Modifications were required to make the algorithm work effectively with the threshold.
- Specifically, when encountering a pivot coefficient equal to 0, it became necessary to set it to η for the modified algorithm to function as expected.

⚖️ **Summary and Recommendations:**
- While LU factorization can be applied to solve linear systems with bounded arithmetic, careful consideration is essential due to numeric instability.
- Tailored measures may be necessary to address specific challenges associated with the algorithm's performance.
  
🚧 **Cautionary Note:**
- The algorithm may not be suitable for scenarios with high condition numbers and infinitesimal coefficients without additional adjustments.

In conclusion, the use of LU factorization presents both opportunities and challenges. Careful parameter tuning and adjustments are crucial to ensure its effectiveness in addressing linear systems with bounded arithmetic. Further research and refinement may be needed to enhance its robustness in the face of diverse scenarios.
