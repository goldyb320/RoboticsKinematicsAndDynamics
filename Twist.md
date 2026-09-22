# What is v? Twists, screws, and the exponential map on SE(3)
 
Study notes for ROB 510 / EECS 567, Lecture 7 (exponential map onto SE(3)), with two interactive visual tools.
 
**Open `index.html` in a browser** to read the notes and use the tools. To host it, enable GitHub Pages for this repo (Settings → Pages → deploy from the main branch, root folder). The page will then be live at `https://goldyb320.github.io/RoboticsKinematicsAndDynamics/`.
 
## What's inside
 
1. The short answer: what v actually is
2. Terms: ω, q, θ, h, v, twist ξ, ξ̂, se(3), SE(3), the exponential map
3. Why a rotating body has no single velocity, plus an interactive velocity-field tool
4. Why v is measured "at the origin," with a worked example
5. Screw motions, plus an interactive 3D screw explorer
6. How it all connects to Props 2.8 and 2.9 in the lecture
7. Why we care: robot joints and product-of-exponentials forward kinematics
8. What v depends on, and common misreadings
9. A quick-reference table
## The one-line version
 
$$v = -\omega \times q + h\,\omega$$
 
v is the velocity of the point at the origin, as if the rigid motion grabbed all of space. ω and q say where the screw axis is, h says how much it slides per turn, and θ says how much it turns. The axis and pitch can be recovered as q = ω × v and h = ωᵀv.
 
## Notes
 
- Everything runs in the browser. Math is rendered with KaTeX, loaded from a CDN, so an internet connection is needed for the equations to display properly.
- The page supports light and dark mode and works on mobile.
