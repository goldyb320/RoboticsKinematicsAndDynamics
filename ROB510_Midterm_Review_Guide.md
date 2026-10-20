# ROB 510 Midterm Review Guide
### Robot Kinematics: rotations → twists → forward kinematics → inverse kinematics → velocities, wrenches, Jacobians

*Built from Lectures 1–14, HW1–HW3, the Winter 2025 and Winter 2026 midterms, and the official midterm review sheet. Exam assumed in person, by hand, no Mathematica.*

---

## 0. How to use this guide (read first, 3 min)

**Exam format (from both past midterms, 70 points, 2 hours):**

| Section | Points | What it is |
|---|---|---|
| Q1–Q5 True/False | 30 (20 statements × 1.5) | One block each: **Rotations**, **Twists**, **Forward Kinematics**, **Inverse Kinematics**, **Velocities/Wrenches/Jacobians** |
| Q6 | 20 | **Forward kinematics** of a drawn robot: twists $\xi_i$, $g_{st}(0)$, D‑H relative twists $\xi_{i-1,i}$, a relative transform, and a twist moved into a new frame with the Adjoint |
| Q7 | 20 | **Velocities, Jacobian, wrench**: spatial Jacobian by hand, the velocity of a point, spatial velocity of a joint motion, joint torques $\tau = J^T F$ |

**What that means for your time.** Half the exam is one long problem on POE twists, and one on Jacobians and velocities. The other 30 points are conceptual T/F questions that hinge on knowing *exactly* what a definition says (for example, "is $\dot R$ in so(3)?", "does Rodrigues' formula need $\|\omega\|=1$?").

**Priority markers used below:**
- ★★★ appears on both exams or is a whole 20‑point problem. Know it cold.
- ★★ appears in T/F or homework. Know the idea and the formula.
- ★ background or proofs. Understand the idea, don't memorize.

**One heads‑up:** the lecture notes stop at the start of wrenches (L14), but both exams test **the manipulator Jacobian** and **$\tau = J^T F$**. That material is probably coming in lectures this week. Section 13 covers it from MLS so you aren't caught off guard.

**Your own practice:** I've deliberately used *different* robots and numbers in the worked examples, so the two old midterms stay fresh for timed practice and HW3 stays yours to do.

```
THE ROADMAP: every idea in the course builds on the one before it

  SO(3) rotations ──(exp of ω̂θ)──▶ Rodrigues' formula
        │                                │
        ▼                                ▼
  SE(3) = (p, R) ──(exp of ξ̂θ)──▶ twists / screws
                                         │
              ┌──────────────────────────┼────────────────────────┐
              ▼                          ▼                        ▼
   Forward kinematics (POE)     Adjoint: moving twists     Inverse kinematics
   g_st = e^ξ1θ1…e^ξnθn g_st(0)  between frames             (Paden–Kahan SP1–3)
              │                          │
              └──────────┬───────────────┘
                         ▼
           Velocities V^s, V^b  →  Jacobian J  →  wrenches, τ = Jᵀ F
```

---

## 1. Rigid body motion and SO(3) (L2–L3) ★★

### 1.1 What "rigid" means

A **rigid body transformation** (RBT) $g:\mathbb R^3\to\mathbb R^3$ must satisfy two things (Def 2.1):

1. **Preserve length:** $\|g(p)-g(q)\| = \|p-q\|$ for all points.
2. **Preserve the cross product:** $g_*(v\times w) = g_*(v)\times g_*(w)$ for all vectors.

Rule 2 is what excludes **reflections**. A mirror preserves lengths but flips handedness, so $x\times y$ would map to $-z$ instead of $z$.

Consequence: inner products are preserved, $v^Tw = g_*(v)^Tg_*(w)$, so orthogonal vectors stay orthogonal.

**Points vs. vectors.** A vector is the difference of two points, $v = q-p$, and $g_*(v) = g(q)-g(p)$. For a general motion $g=(p,R)$, a vector only gets **rotated**: $g_*(v) = Rv$. Translation cancels out in the difference. This matters later for homogeneous coordinates.

### 1.2 Rotation matrices: what the columns mean ★★★

$$R_{ab} = \begin{bmatrix} x_{ab} & y_{ab} & z_{ab}\end{bmatrix}$$

**Column $i$ of $R_{ab}$ is B's $i$‑th axis written in A's coordinates.** This one fact lets you write any rotation by looking at a picture.

```
        z_A                       Example: B is A rotated +90° about z_A
         │                         x_B points along +y_A  → column 1 = (0,1,0)
         │      y_B                y_B points along −x_A  → column 2 = (−1,0,0)
         │     ╱                   z_B = z_A              → column 3 = (0,0,1)
         │    ╱
         └──────── y_A                       ┌ 0 −1  0 ┐
        ╱  ╲                         R_ab =  │ 1  0  0 │ = R_z(90°)
       ╱    ╲ x_B (= +y_A dir)               └ 0  0  1 ┘
     x_A
```

**What $R_{ab}$ does, two equivalent readings:**
- **Change of coordinates:** $q_a = R_{ab}\,q_b$. Give it a point in B's coordinates and get the same point in A's coordinates.
- **Active rotation:** it rotates a vector, relative to a fixed frame.

**Composition and inverse (subscripts "cancel" like fractions):**
$$R_{ac} = R_{ab}R_{bc},\qquad R_{ba} = R_{ab}^{-1} = R_{ab}^T$$

> **T/F pattern ★★★:** "A has orientation $R_1$ relative to B, and C has orientation $R_2$ relative to B. What is A relative to C?" Translate the words into subscripts: $R_1 = R_{BA}$ and $R_2 = R_{BC}$. Then $R_{CA} = R_{CB}R_{BA} = R_2^TR_1$. Always convert words into subscripts and chain them so the inner letters match.

### 1.3 SO(3), groups, and abelian groups ★★

$$SO(3) = \{R\in\mathbb R^{3\times3} : RR^T = R^TR = I,\ \det R = +1\}$$

- $RR^T=I$ (orthogonal) gives $\det R = \pm1$.
- "Special" means $\det = +1$, which excludes reflections.
- $\det R_{ab} = x_{ab}^T(y_{ab}\times z_{ab}) = x^Tx = 1$, because the frame is right‑handed.

**Group axioms** (a set with an operation $\circ$):
1. Closure
2. Identity
3. Inverse
4. Associativity

**SO(3) is a group under matrix multiplication** (proof in L2):
- closure: $(R_1R_2)^T(R_1R_2) = I$ and $\det = 1\cdot1$
- identity: $I$
- inverse: $R^T$
- associativity: inherited from matrix multiplication

**Abelian** means every pair commutes: $g_1g_2 = g_2g_1$.

| Set | Operation | Group? | Abelian? | Why |
|---|---|---|---|---|
| $SO(2)$ | multiplication | yes | **yes** | Every element is $R(\theta)$, and $R(a)R(b) = R(a+b) = R(b)R(a)$ |
| $SO(3)$, $SO(n\ge3)$ | multiplication | yes | **no** | 3D rotations about different axes don't commute (order matters) |
| $so(3)$ (skew matrices) | **addition** | yes | **yes** | It's a vector space, and matrix addition always commutes |
| $SE(3)$ | multiplication | yes | no | contains SO(3) |
| $(\mathbb Z,+)$, $(\mathbb R,+)$ | addition | yes | yes | |

> **Trap:** so(3) is **not** closed under multiplication (the product of two skew matrices is usually not skew), so it is not a group under multiplication. Under addition it is an abelian group. Read which operation the question names.

### 1.4 The hat map and cross product ★★★

$$a\times b = \hat a\, b,\qquad \hat a = \begin{bmatrix}0&-a_3&a_2\\ a_3&0&-a_1\\ -a_2&a_1&0\end{bmatrix},\qquad \hat a^T = -\hat a$$

**Memory trick:** zeros on the diagonal; $+a_1$ sits at position (3,2), $+a_2$ at (1,3), $+a_3$ at (2,1), and each mirror entry is the negative.

**Vee map:** $(\hat a)^\vee = a$. You read $a$ back off as $(\hat a_{32}, \hat a_{13}, \hat a_{21})$.

**Key facts:**
- $\widehat{Rw} = R\,\hat w\,R^T$. Rotating the axis is the same as conjugating the skew matrix.
- $\hat w\,w = w\times w = 0$, so **$w$ is in the null space of $\hat w$**.
- **$\hat w$ is never invertible.** Its determinant is always 0, because $\hat w w = 0$ and its eigenvalues are $0, \pm i\|w\|$.

> **Trap ★★:** "If $\|\omega\|=1$, then $\hat\omega$ is invertible with $\det\hat\omega = 1$." **False.** A 3×3 skew matrix always has $\det = 0$. (The determinant of an odd‑dimension skew matrix is $\det(\hat w) = \det(\hat w^T) = \det(-\hat w) = (-1)^3\det\hat w$, so it must be 0.)

---

## 2. Euler angles (L3–L4) ★★

### 2.1 Elementary rotations (memorize these)

$$R_x(\phi)=\begin{bmatrix}1&0&0\\0&c&-s\\0&s&c\end{bmatrix}\quad
R_y(\beta)=\begin{bmatrix}c&0&s\\0&1&0\\-s&0&c\end{bmatrix}\quad
R_z(\alpha)=\begin{bmatrix}c&-s&0\\s&c&0\\0&0&1\end{bmatrix}$$

**Pattern:** the axis you rotate about gets the "1" row and column. The other 2×2 block is a 2D rotation with $-s$ in the upper right, **except $R_y$, where the $-s$ is in the lower left**. That's because the cyclic order is $z\to x$ for a y‑rotation.

### 2.2 ZYZ and ZYX (yaw–pitch–roll)

**Intrinsic** (rotate about the body's *new* axes each time) means you **post‑multiply** in the order you say them:

- **ZYZ:** rotate about z by α, then the new y by β, then the newest z by γ:
$$R_{ab} = R_z(\alpha)R_y(\beta)R_z(\gamma)=\begin{bmatrix}c_\alpha c_\beta c_\gamma - s_\alpha s_\gamma & -c_\alpha c_\beta s_\gamma - s_\alpha c_\gamma & c_\alpha s_\beta\\ s_\alpha c_\beta c_\gamma + c_\alpha s_\gamma & -s_\alpha c_\beta s_\gamma + c_\alpha c_\gamma & s_\alpha s_\beta\\ -s_\beta c_\gamma & s_\beta s_\gamma & c_\beta\end{bmatrix}$$
- **ZYX (yaw φ, pitch θ, roll ψ):** $R = R_z(\phi)R_y(\theta)R_x(\psi)$. Its bottom row is $(-s_\theta,\ c_\theta s_\psi,\ c_\theta c_\psi)$ and its first column is $(c_\phi c_\theta, s_\phi c_\theta, -s_\theta)$.

**Extrinsic** (about the fixed world axes) means you **pre‑multiply**. The same three rotations done about fixed axes give the product in reverse order. The intrinsic sequence z, y′, z″ equals the extrinsic sequence z, y, z with the matrices multiplied the other way around.

### 2.3 Inverse: from R back to angles

**ZYZ**, when $\sin\beta\neq0$:
$$\beta = \operatorname{atan2}\!\big(\sqrt{r_{31}^2+r_{32}^2},\ r_{33}\big),\quad
\alpha = \operatorname{atan2}(r_{23}/s_\beta,\ r_{13}/s_\beta),\quad
\gamma = \operatorname{atan2}(r_{32}/s_\beta,\ -r_{31}/s_\beta)$$

**Where these come from:** look at the third column, $(c_\alpha s_\beta, s_\alpha s_\beta, c_\beta)$, and the third row, $(-s_\beta c_\gamma, s_\beta s_\gamma, c_\beta)$. You don't need to memorize the formulas, just how to read them off the matrix.

**ZYX:** $\theta = \operatorname{atan2}(-r_{31}, \sqrt{r_{11}^2+r_{21}^2})$, $\phi = \operatorname{atan2}(r_{21}, r_{11})$, $\psi = \operatorname{atan2}(r_{32}, r_{33})$.

### 2.4 Singularities and the topology fact ★★★

- **ZYZ is singular at $\beta = 0$ or $\pi$.** At $\beta=0$ the first and last rotations are about the *same* axis, so only $\alpha+\gamma$ matters. For $R = I$, every triple $(\alpha, 0, -\alpha)$ works, which gives **infinitely many solutions**.
- **ZYX is singular at pitch $\theta = \pm\pi/2$** (the plane pointing straight down). This is gimbal lock.
- **Topological fact:** *no* 3‑parameter representation of SO(3) is free of singularities. The 9‑element matrix is global but redundant. Quaternions (4 numbers) are global.

```
   Gimbal lock (ZYZ, β = 0):

   ═══[ z‑axis ]═══(  y  )═══[ z‑axis ]═══
       rotation 1     β=0     rotation 3      ← axes 1 and 3 line up,
                                                 so α and γ do the same thing
```

> **T/F ★★★:** "The map from Euler angles $(\alpha,\beta,\gamma)$ to $R\in SO(3)$ is injective." **False.** Many triples give the same R: the singularity at β=0, plus 2π periodicity, plus the $(\alpha+\pi, -\beta, \gamma+\pi)$ twin.
> Euler angles are a **local** parameterization, not a global one.

---

## 3. Exponential coordinates of rotation (L4–L5, HW1) ★★★

### 3.1 Where the exponential comes from

Spin a point $q$ about a fixed unit axis $\omega$ through the origin at unit speed. Its velocity is
$$\dot q(t) = \omega\times q(t) = \hat\omega\,q(t).$$
This is a linear time‑invariant ODE, so its solution is $q(t) = e^{\hat\omega t}q(0)$. After rotating by angle $\theta$:
$$\boxed{R(\omega,\theta) = e^{\hat\omega\theta}}$$

```
        ω (unit axis)
        ↑
        │      q̇ = ω × q   (tangent to the circle)
        │    ↗
        │  ● q(t)
        │ ╱
     ───┼──────  the point sweeps a circle around ω;
        │          the angle swept is θ
```

**$\omega\theta\in\mathbb R^3$ are the exponential coordinates of R.** That's three numbers: the axis direction scaled by the angle.

### 3.2 so(3) and Lemma 2.3 (powers of $\hat a$) ★★

$so(3) = \{S\in\mathbb R^{3\times3}: S^T = -S\}$, a 3‑dimensional real vector space.

**Lemma 2.3:**
$$\hat a^2 = aa^T - \|a\|^2I,\qquad \hat a^3 = -\|a\|^2\hat a$$

With $a = \omega\theta$ and $\|\omega\|=1$, the powers cycle:

| power | value |
|---|---|
| $\hat a$ | $\theta\hat\omega$ |
| $\hat a^2$ | $\theta^2\hat\omega^2$ |
| $\hat a^3$ | $-\theta^3\hat\omega$ |
| $\hat a^4$ | $-\theta^4\hat\omega^2$ |
| $\hat a^5$ | $+\theta^5\hat\omega$ |

So the odd powers are $(-1)^{n-1}\theta^{2n-1}\hat\omega$ and the even powers are $(-1)^{n+1}\theta^{2n}\hat\omega^2$. These are exactly HW1's Lemmas 0.1 and 0.2, proved by induction.

### 3.3 Rodrigues' formula ★★★

Plug the cycle into $e^{\hat\omega\theta} = I + \hat\omega\theta + \frac{(\hat\omega\theta)^2}{2!}+\cdots$ and group the terms:
$$e^{\hat\omega\theta} = I + \underbrace{\Big(\theta - \tfrac{\theta^3}{3!} + \cdots\Big)}_{\sin\theta}\hat\omega + \underbrace{\Big(\tfrac{\theta^2}{2!} - \tfrac{\theta^4}{4!} + \cdots\Big)}_{1-\cos\theta}\hat\omega^2$$
$$\boxed{e^{\hat\omega\theta} = I + \hat\omega\sin\theta + \hat\omega^2(1-\cos\theta)}\qquad(\|\omega\| = 1\ \text{required!})$$

> **T/F ★★★:** "Rodrigues' formula does not hold for some $\hat\omega\in so(3)$ and $\theta$." Rodrigues' formula **only holds when $\|\omega\|=1$**. The derivation used $\hat\omega^3 = -\hat\omega$, which needs a unit vector. For a non‑unit $\omega$, first normalize: $\omega\theta = \frac{\omega}{\|\omega\|}(\|\omega\|\theta)$.

**By‑hand tip:** for a coordinate axis you never need Rodrigues. $e^{\hat z\theta} = R_z(\theta)$, $e^{\hat x\theta} = R_x(\theta)$, $e^{\hat y\theta} = R_y(\theta)$. Use Rodrigues only for tilted axes. Also useful: $\hat\omega^2 = \omega\omega^T - I$.

### 3.4 exp lands in SO(3) (Prop 2.4) ★

- **Orthogonal:** $(e^{\hat\omega\theta})^{-1} = e^{-\hat\omega\theta} = e^{\hat\omega^T\theta} = (e^{\hat\omega\theta})^T$.
- **$\det = +1$:** orthogonality gives $\det = \pm1$. At $\theta=0$, $\det I = 1$, and $\det$ is continuous in $\theta$, so it can't jump to $-1$.

### 3.5 Surjectivity: every R has an axis‑angle (Prop 2.5, Euler's Thm 2.6) ★★★

Write out Rodrigues' formula entrywise, with $v_\theta = 1-\cos\theta$:
$$R = \begin{bmatrix}\omega_1^2v_\theta + c_\theta & \omega_1\omega_2v_\theta - \omega_3s_\theta & \omega_1\omega_3v_\theta+\omega_2s_\theta\\ \omega_1\omega_2v_\theta+\omega_3s_\theta & \omega_2^2v_\theta+c_\theta & \omega_2\omega_3v_\theta - \omega_1s_\theta\\ \omega_1\omega_3v_\theta - \omega_2s_\theta & \omega_2\omega_3v_\theta+\omega_1s_\theta & \omega_3^2v_\theta+c_\theta\end{bmatrix}$$

**Step 1, the angle from the trace:**
$$\operatorname{tr}R = 1+2\cos\theta\ \Rightarrow\ \boxed{\theta = \cos^{-1}\!\Big(\frac{\operatorname{tr}R-1}{2}\Big)}$$

**Step 2, the axis from the skew part:**
$$\boxed{\omega = \frac{1}{2\sin\theta}\begin{bmatrix}r_{32}-r_{23}\\ r_{13}-r_{31}\\ r_{21}-r_{12}\end{bmatrix}}\quad(\sin\theta\neq0)$$

**The two special cases (HW1 P2, P3):**

| Case | What happens | How to get ω |
|---|---|---|
| $\theta = 0$ ($R = I$) | A **true singularity**. Any ω works, since $I = e^{\hat\omega\cdot0}$. The axis is meaningless when you haven't rotated. | arbitrary |
| $\theta = \pi$ | **Not** a singularity. The formula breaks (divide by $\sin\pi = 0$), but ω is still defined, up to sign. | $R = I + 2\hat\omega^2 = 2\omega\omega^T - I$, so $\omega\omega^T = \frac{R+I}{2}$. Then $\omega_i = \pm\sqrt{\frac{r_{ii}+1}{2}}$, with the relative signs taken from $r_{ij} = 2\omega_i\omega_j$. |

**Not injective.** These all give the same R:
- $(\omega,\theta)$
- $(-\omega, -\theta)$
- $(-\omega, 2\pi-\theta)$
- $(\omega, \theta+2\pi k)$

> **T/F ★★★:** "There exists R with infinitely many exponential coordinates $(\omega,\theta)$, $\theta\in[0,2\pi)$." **True:** $R = I$, where $\theta = 0$ and any ω works.
> **"exp: so(3) → SO(3) is surjective"** is true. **"Injective"** is false.

### 3.6 Eigenvalues (HW1 P4) ★★

- $\hat\omega$ (unit ω) has eigenvalues $0, \pm i$. The eigenvector for 0 is ω itself, since $\hat\omega\omega = 0$.
- $R = e^{\hat\omega\theta}$ has eigenvalues $e^{0} = 1, e^{\pm i\theta}$, **with the same eigenvectors**. If $\hat\omega x = \lambda x$, then $e^{\hat\omega\theta}x = e^{\lambda\theta}x$, term by term in the series.
- **Geometric meaning:** the eigenvector with eigenvalue 1 is the **rotation axis**. It's the one direction R leaves alone.

> **Classic T/F ★★★:** "$A^{-1} = A^T$ and the eigenvalues are $1, e^{\pm i\pi/2}$. Is $A\in SO(3)$?"
> Orthogonality gives $\det = \pm1$, and $\det$ = the product of the eigenvalues $= 1\cdot e^{i\pi/2}e^{-i\pi/2} = 1$. So **yes**.
> General rule: **an orthogonal matrix is in SO(3) exactly when the product of its eigenvalues is +1.** A reflection would have an eigenvalue of −1 where you'd expect the 1.

### 3.7 SO(2) (HW1 P5–P6) ★

- $SO(2) = \{R(\theta)\}\cong S^1$, the circle: one parameter with period 2π.
- $\hat\omega = \begin{bmatrix}0&-\omega\\\omega&0\end{bmatrix}$ satisfies $\hat\omega^2 = -\omega^2I$, so $e^{\hat\omega\theta} = \cos(\omega\theta)I + \sin(\omega\theta)\hat\omega/\omega = R(\omega\theta)$.
- **exp: so(2)→SO(2) is surjective but not injective** ($\omega\theta$ and $\omega\theta+2\pi$ give the same R).
- **SO(2) commutes:** $R\hat\omega R^T = \hat\omega$, $Re^{\hat\omega\theta}R^T = e^{\hat\omega\theta}$, and $\frac{d}{dt}e^{\hat\omega\theta} = \hat\omega\dot\theta e^{\hat\omega\theta} = e^{\hat\omega\theta}\hat\omega\dot\theta$.
- **In 3D this fails.** $R\hat\omega R^T = \widehat{R\omega}$, which is a *different* axis.

---

## 4. SE(3) and homogeneous coordinates (L5–L6) ★★★

### 4.1 Configuration of a body

The configuration of frame B relative to A is $g_{ab} = (p_{ab}, R_{ab})$:
- $p_{ab}$ is the position of B's origin, in A's coordinates.
- $R_{ab}$ is B's axes, in A's coordinates.

$$SE(3) = \mathbb R^3\times SO(3),\qquad q_a = p_{ab} + R_{ab}q_b$$

**Elements of SE(3) are RBTs (Prop 2.7):**
- Lengths: $\|g(p_1) - g(p_2)\| = \|R(p_1-p_2)\|$.
- Cross products: $g_*(v) = Rv$ for vectors.

### 4.2 Homogeneous form (this is what you'll compute with)

$$\bar g_{ab} = \begin{bmatrix}R_{ab}&p_{ab}\\0&1\end{bmatrix}\in\mathbb R^{4\times4},\qquad \bar q = \begin{bmatrix}q\\1\end{bmatrix}\ (\text{point}),\qquad \bar v = \begin{bmatrix}v\\0\end{bmatrix}\ (\text{vector})$$

**The bottom entry is the type check:**

| Expression | Bottom entry | Result |
|---|---|---|
| vector ± vector | 0 | vector |
| point + vector | 1 | point |
| point − point | 0 | vector |
| point + point | 2 | **meaningless** |

A vector ignores the translation: $\bar g\bar v = (Rv, 0)$.

**Composition and inverse:**
$$\bar g_{ac} = \bar g_{ab}\bar g_{bc} = \begin{bmatrix}R_{ab}R_{bc} & R_{ab}p_{bc}+p_{ab}\\0&1\end{bmatrix},\qquad \boxed{\bar g^{-1} = \begin{bmatrix}R^T & -R^Tp\\0&1\end{bmatrix}}$$

**$\det\bar g = \det R = 1$, always**, even for pure translation. This is block upper triangular, so the determinant is $\det R\cdot1$.

> **T/F ★★:** "For pure translation g, $\det g = 0$." **False.** It's 1. (Don't confuse g with the twist $\hat\xi$, which *is* singular, with $\det\hat\xi = 0$.)

---

## 5. Twists and the exponential map on SE(3) (L6–L7) ★★★

### 5.1 Twist for a rotation about an *offset* axis

Let ω be a unit axis and **q any point on the axis**. A point $p$ spinning around that axis moves with velocity
$$\dot p = \omega\times(p-q) = \hat\omega p - \hat\omega q = \hat\omega p + v,\qquad v := -\omega\times q = q\times\omega$$

In homogeneous form this is linear: $\dot{\bar p} = \hat\xi\,\bar p$, so $\bar p(\theta) = e^{\hat\xi\theta}\bar p(0)$, where
$$\hat\xi = \begin{bmatrix}\hat\omega & v\\0&0\end{bmatrix}\in se(3)\ (4\times4),\qquad \xi = \begin{bmatrix}v\\\omega\end{bmatrix}\in\mathbb R^6\ \text{(twist coordinates; } v \text{ on top in MLS)}$$

```
            ω                     v = −ω × q
            ↑                     • depends on where the axis IS (through q)
            │ ↻                   • does NOT depend on which q you pick on the axis:
            │                       q + tω gives −ω×(q + tω) = −ω×q, the same v
    q ●─────┤
            │                     • axis through the origin (q = 0) gives v = 0
            │
```

**Bottom row of $\hat\xi$ is all zeros** (including the bottom‑right entry). Compare $\bar g$, whose bottom‑right entry is 1. The twist is a *velocity‑like* object, and a velocity times a point gives a vector:
$$\hat\xi\bar p = \begin{bmatrix}\omega\times p + v\\0\end{bmatrix}$$
That bottom 0 marks it as a vector.

### 5.2 Twist types ★★★ (know these cold)

| Joint / motion | Twist coordinates $\xi = (v,\omega)$ | Notes |
|---|---|---|
| **Revolute** (pure rotation, pitch 0) | $\begin{bmatrix}-\omega\times q\\\omega\end{bmatrix}$, $\|\omega\| = 1$, q on the axis | the axis through the origin gives $v=0$ |
| **Prismatic** (pure translation, pitch ∞) | $\begin{bmatrix}v\\0\end{bmatrix}$, $\|v\| = 1$ | **where** the slider sits doesn't matter, only its direction |
| **Screw/helical** (pitch h) | $\begin{bmatrix}-\omega\times q + h\omega\\\omega\end{bmatrix}$ | translate $h\theta$ along ω while rotating θ |

**Fast cross products you'll use constantly** (right‑hand cyclic $x\to y\to z\to x$):
- $z\times x = y$, $x\times y = z$, $y\times z = x$
- reversed order flips the sign: $x\times z = -y$, and so on
- $-\omega\times q = q\times\omega$, often easier to remember

> **Worked mini‑example.** Revolute about +x through $q = (0, a, b)$:
> $v = q\times\omega = (0,a,b)\times(1,0,0) = (a\cdot0 - b\cdot0,\ b\cdot1 - 0\cdot0,\ 0\cdot0 - a\cdot1) = (0, b, -a)$.
> So $\xi = (0, b, -a, 1, 0, 0)$.
> Notice the x‑coordinate of q never matters (it's along the axis), as expected.

> **T/F ★★:** "A unit twist has 2‑norm 1." **False** in general. *Unit* means $\|\omega\| = 1$ (or $\|v\|=1$ if $\omega=0$). The full vector has $\|\xi\| = \sqrt{\|v\|^2+1}$, which is 1 only when the axis passes through the origin.

### 5.3 The exponential map formulas ★★★

**Pure translation ($\omega = 0$).** Since $\hat\xi^2 = 0$, the series stops:
$$e^{\hat\xi\theta} = \begin{bmatrix}I & v\theta\\0&1\end{bmatrix}$$

**Rotation or screw ($\omega\neq0$, $\|\omega\|=1$):**
$$\boxed{e^{\hat\xi\theta} = \begin{bmatrix}e^{\hat\omega\theta} & (I - e^{\hat\omega\theta})(\omega\times v) + \omega\omega^Tv\,\theta\\0&1\end{bmatrix}}$$

**The by‑hand shortcut for revolute joints** (the one you actually want on the exam):
$$\boxed{\text{pure rotation about an axis through } q:\quad e^{\hat\xi\theta} = \begin{bmatrix}e^{\hat\omega\theta} & (I - e^{\hat\omega\theta})\,q\\0&1\end{bmatrix}}$$

**Why the shortcut is the same thing:** with $v = -\omega\times q$, we get $\omega\times v = q - (\omega^Tq)\omega = q_\perp$, the part of q perpendicular to the axis. Also $\omega^Tv = 0$. And $(I-R)q = (I-R)q_\perp$, because R leaves the along‑axis part alone. Pick the **simplest q on the axis** and you're done in one line.

**Geometric reading:** rotate by R about the origin, then shift by $(I-R)q$ to put the pivot back at q:
$$p \mapsto R(p - q) + q$$

**Screw ($h\neq0$):** add $h\theta\,\omega$ to the top‑right entry:
$$e^{\hat\xi\theta} = \begin{bmatrix}R & (I - R)q + h\theta\omega\\0&1\end{bmatrix}$$

> **T/F ★★★:** "The general formula $e^{\hat\xi\theta} = \begin{bmatrix}e^{\hat\omega\theta}&(I-e^{\hat\omega\theta})(\omega\times v)+\omega\omega^Tv\theta\\0&1\end{bmatrix}$ works for any twist, including pure translation." **False.** Set $\omega = 0$ and the top‑right entry is $0$, but the correct answer is $v\theta$. The formula also assumed $\|\omega\| = 1$. Translation needs its own case.

> **Worked example (L10, SCARA joint 2).** $\omega = z$, $q = (0, l_1, 0)$:
> $$e^{\hat\xi_2\theta_2} = \begin{bmatrix}R_z(\theta_2) & (I - R_z(\theta_2))\begin{bmatrix}0\\l_1\\0\end{bmatrix}\\0&1\end{bmatrix} = \begin{bmatrix}c_2&-s_2&0& l_1s_2\\ s_2&c_2&0& l_1(1-c_2)\\0&0&1&0\\0&0&0&1\end{bmatrix}$$
> **Sanity check:** a point on the axis, $(0, l_1, 0)$, must stay put. $R_z(0,l_1,0) + (l_1s_2, l_1-l_1c_2, 0) = (-l_1s_2 + l_1s_2,\ l_1c_2 + l_1 - l_1c_2,\ 0) = (0, l_1, 0)$ ✓.
> **Always do this check.** It catches sign errors in seconds.

### 5.4 How the formula is proven (L7, the similarity trick) ★

1. Move the frame to the closest point on the axis using $g = \begin{bmatrix}I & \omega\times v\\0&1\end{bmatrix}$.
2. In the new coordinates, the twist becomes $\hat\xi' = g^{-1}\hat\xi g = \begin{bmatrix}\hat\omega & h\omega\\0&0\end{bmatrix}$, an axis through the origin. Its exponential is easy: $\begin{bmatrix}e^{\hat\omega\theta} & h\omega\theta\\0&1\end{bmatrix}$.
3. Transform back: $e^{\hat\xi\theta} = g\,e^{\hat\xi'\theta}g^{-1}$.

The idea worth keeping: **$e^{g\hat\xi g^{-1}\theta} = g\,e^{\hat\xi\theta}g^{-1}$.** This identity is the whole reason the Adjoint works (Section 8).

### 5.5 Surjectivity onto SE(3) (Prop 2.9) ★★

**Given any $g = (R, p)$, there is a twist and θ with $g = e^{\hat\xi\theta}$:**
- **$R = I$** (pure translation): $v = p/\|p\|$, $\theta = \|p\|$.
- **$R\neq I$:** first solve $e^{\hat\omega\theta} = R$ (Prop 2.5). Then solve $Av = p$ with $A = (I - e^{\hat\omega\theta})\hat\omega + \omega\omega^T\theta$. A is invertible for $\theta\neq0$, because its two pieces have mutually orthogonal null spaces.

**exp: se(3) → SE(3) is surjective, not injective** (inherited from SO(3)).
**Exponential coordinates of g:** $\xi\theta\in\mathbb R^6$.

### 5.6 Planar version, se(2) (from the review sheet) ★

$$\hat\xi = \begin{bmatrix}\hat\omega & v\\0&0\end{bmatrix}\in\mathbb R^{3\times3},\quad \hat\omega = \begin{bmatrix}0&-\omega\\\omega&0\end{bmatrix}\in so(2),\quad v\in\mathbb R^2,\quad \xi = (v,\omega)\in\mathbb R^3$$

A rotation about the planar point $q = (q_x,q_y)$ has $v = -\hat\omega q = (\omega q_y, -\omega q_x)$.

**Every planar rigid motion is a pure rotation** about some point (possibly infinitely far away, which is the pure‑translation case).

---

## 6. Screws (L8) ★★

**Q:** for an arbitrary twist $(v,\omega)$, can we always find q with $v = -\omega\times q$?
**A: No.** $-\omega\times q$ is always ⟂ ω. If v has a component along ω, you need a **pitch** term:
$$v = -\omega\times q + h\omega,\qquad \boxed{h = \frac{\omega^Tv}{\|\omega\|^2}}\ \ (\text{translation per radian})$$

```
   Screw motion = rotate θ about the axis  +  slide hθ along the same axis

        q ●━━━━━━━━━━━━━━━━━━━━━━━━━━▶ ω  (screw axis)
             ↻ θ          ──hθ──▶

   h = 0 : pure rotation      h = ∞ : pure translation (ω = 0)
```

**Screw ↔ twist conversion:**

| Screw attribute | From the twist $\xi = (v,\omega)$ |
|---|---|
| Pitch | $h = \omega^Tv/\|\omega\|^2$ ($\infty$ if $\omega=0$) |
| Axis | $\{\frac{\omega\times v}{\|\omega\|^2} + \lambda\omega\}$ if $\omega\neq0$; $\{0 + \lambda v\}$ if $\omega=0$ |
| Magnitude | $M = \|\omega\|$ if $\omega\neq0$, else $\|v\|$ |

The point on the axis closest to the origin is $\omega\times v/\|\omega\|^2$.

**Unit twist and magnitude.** If $\xi_{unit}$ has $\|\omega\| = 1$ (or $\|v\| = 1$), then for $\xi = \xi_{unit}\theta$ the magnitude is $M = \|\omega\theta\| = |\theta|$. So **θ plays the role of the screw magnitude** (its absolute value, if you want to be picky about signs).

**Chasles (Thm 2.11):** *every* rigid motion in 3D is a screw motion. Instantaneous velocities are twists too, so **every velocity $\hat V^s$ or $\hat V^b$ has a screw** with zero, finite, or infinite pitch.

> **Worked example (HW2 P1 style, a different bolt).** A nut threads onto a bolt parallel to z through $(2, 0, 0)$ with 8 threads per inch.
> - One turn ($2\pi$ rad) moves it $1/8$ inch, so $h = \frac{1}{8\cdot2\pi} = \frac{1}{16\pi}$ in/rad.
> - Sign: for a right‑hand thread turned the right‑hand way about $+z$, it advances along $+z$.
>
> $$\xi = \begin{bmatrix}-z\times q + hz\\ z\end{bmatrix},\quad -z\times(2,0,0) = -(0,2,0)\ \Rightarrow\ \xi = \Big(0,\ -2,\ \tfrac{1}{16\pi},\ 0,0,1\Big)$$

---

## 7. Forward kinematics with the Product of Exponentials (L8–L10) ★★★ (Q6)

### 7.1 The setup

- **Joint space** $Q$: one entry per joint ($S^1$ for a revolute joint, $\mathbb R$ for a prismatic one). Example: SCARA is $S^1\times S^1\times S^1\times\mathbb R = T^3\times\mathbb R$.
- **Joint i connects link i−1 to link i.** S is attached to the base (link 0), and T to the last link.
- **The forward kinematics map** is $g_{st}: Q\to SE(3)$.

### 7.2 The four‑step POE recipe ★★★

1. **Reference configuration** $\theta = 0$: write $g_{st}(0) = \begin{bmatrix}R_{st}(0) & p_{st}(0)\\0&1\end{bmatrix}$.
   - $p$ is the tool origin measured from S.
   - The columns of R are the tool axes in S coordinates.
   - Frames drawn parallel give $R = I$.
2. **Twist for each joint**, in S coordinates, at $\theta = 0$:
   - revolute: $\xi_i = (-\omega_i\times q_i, \omega_i)$
   - prismatic: $\xi_i = (v_i, 0)$
3. **Exponentials** $e^{\hat\xi_i\theta_i}$, using the shortcut from 5.3.
4. **Multiply:**
$$\boxed{g_{st}(\theta) = e^{\hat\xi_1\theta_1}e^{\hat\xi_2\theta_2}\cdots e^{\hat\xi_n\theta_n}\,g_{st}(0)}$$

**The same recipe works for any link frame or point:**
- $g_{sl_i}(\theta) = e^{\hat\xi_1\theta_1}\cdots e^{\hat\xi_i\theta_i}g_{sl_i}(0)$, using **only the joints before that link**.
- A point on link i: $\bar p(\theta) = e^{\hat\xi_1\theta_1}\cdots e^{\hat\xi_i\theta_i}\bar p(0)$.

### 7.3 Reading a figure: the procedure that gets you points ★★★

```
For each joint i:
  ┌──────────────────────────────────────────────────────────────────┐
  │ 1. Revolute or prismatic? (cylinder with a rotation arrow vs.    │
  │    box/slider with a straight arrow)                             │
  │ 2. Direction ω_i (or v_i): read the arrow and use the RIGHT-HAND │
  │    RULE. Thumb along ω and fingers curl in the +θ direction.     │
  │    If the drawn arrow curls the "wrong" way, ω points the other  │
  │    way (−x instead of +x). The exam says "use the angle          │
  │    conventions shown," so obey the arrows.                       │
  │ 3. A point q_i on the axis: choose the EASIEST point (the one    │
  │    with the most zeros). Lengths come from the dimension lines.  │
  │ 4. ξ_i = (q_i × ω_i, ω_i)   or   (v_i, 0)                        │
  │ 5. Sanity: is v_i ⟂ ω_i? (it must be for revolute)               │
  └──────────────────────────────────────────────────────────────────┘
```

**Axis conventions in MLS figures:** z up, y to the right, **x out of the page**. In a drawing, an axis that goes diagonally "into the page" (up and to the right) is **−x**. One that comes diagonally out (down and to the left) is **+x**.

> **Worked example (4‑DOF, my own robot; draw it as you read).**
>
> ```
>     z                             J2: revolute about +x at (0,0,l0)
>     ↑      J2 ⊙──── J3 ═══▶ J4 ⊙──────── T     J3: prismatic along +y
>     │      (x)      (slide y)  (x)             J4: revolute about +x at (0,l1,l0)
>     │ J1 (about z)   ←─ l1 ─→   ←── l2 ──→     T at (0, l1+l2, l0), frames parallel to S
>     │  │l0
>     S────────▶ y
> ```
>
> - $g_{st}(0) = \begin{bmatrix}I & (0,\ l_1+l_2,\ l_0)^T\\0&1\end{bmatrix}$
> - $\xi_1$: $\omega = z$, $q = 0$, so $\xi_1 = (0,0,0,\ 0,0,1)$.
> - $\xi_2$: $\omega = x$, $q = (0,0,l_0)$. $v = q\times x = (0,0,l_0)\times(1,0,0) = (0,\ l_0,\ 0)$, so $\xi_2 = (0,l_0,0,\ 1,0,0)$.
> - $\xi_3$: prismatic, so $\xi_3 = (0,1,0,\ 0,0,0)$.
> - $\xi_4$: $\omega = x$, $q = (0,l_1,l_0)$. $v = (0,l_1,l_0)\times(1,0,0) = (0,\ l_0,\ -l_1)$, so $\xi_4 = (0,l_0,-l_1,\ 1,0,0)$.
>
> **Check:** $v_4\cdot\omega_4 = 0$ ✓. Moving q along x changes nothing ✓.

**Small robots you should be able to write instantly** (from lecture and review):
- **SCARA** (L10): $\xi_{1,2,3} = (q_i\times z, z)$ with $q_1 = 0$, $q_2 = (0,l_1,0)$, $q_3 = (0, l_1+l_2, 0)$. This gives $\xi_2 = (l_1,0,0,0,0,1)$, $\xi_3 = (l_1+l_2,0,0,0,0,1)$, and $\xi_4 = (0,0,1,0,0,0)$. Also $g_{st}(0) = [I \mid (0, l_1+l_2, l_0)]$.
- **Elbow manipulator** (MLS Fig 3.11): $\omega_1 = z$; $\omega_2 = \omega_3 = -x$; $\omega_4 = z$; $\omega_5 = -x$; $\omega_6 = y$. Points: $q_1 = q_2 = (0,0,l_0)$, $q_3 = (0,l_1,l_0)$, $q_{4,5,6} = (0, l_1+l_2, l_0)$. The wrist axes meet at one point.

### 7.4 Order and physical motion ★★★ (T/F favorite)

**Order of exponentials:**
- The POE formula reads as "move the **last** joint first, then the next to last, …, the first joint last." Each $e^{\hat\xi_i\theta_i}$ uses the twist **at the reference configuration**. That's valid because when you move joint n first, the axes of joints 1…n−1 haven't moved yet.
- **Physically, order doesn't matter.** The final pose depends only on the joint values, not on the order the motors actually moved. The math just has to account for it consistently (Section 8.3).

**Exponentials don't add:**
$$e^{\hat\xi_1}e^{\hat\xi_2}e^{\hat\xi_3}\neq e^{\hat\xi_1+\hat\xi_2+\hat\xi_3}$$
in general. $e^Ae^B = e^{A+B}$ only if A and B commute, and twists usually don't. (SE(3) is not abelian.)

**Rigid motions preserve distance.** $\|g\bar x\| = \|x\|$ for a vector x, and the difference of two points is a vector. So for example:
$$\big\|e^{\hat\xi_1\theta_1}e^{\hat\xi_2\theta_2}\big(e^{\hat\xi_3\theta_3}\bar p - \bar q\big)\big\| = \big\|e^{\hat\xi_3\theta_3}\bar p - \bar q\big\|$$
This is exactly the trick that drives SP3 in inverse kinematics.

**$e^{\hat\xi\theta}g = g$ is generally false.** Exponentials don't leave frames alone. They fix only points on their own axis (for a revolute joint).

**A tool point vs. a tool frame:**
- A point $p_t$ given in **tool coordinates** maps to its current spatial position by $g_{st}(\theta)\bar p_t$.
- A point given in **spatial coordinates at the reference configuration** moves by $e^{\hat\xi_1\theta_1}\cdots e^{\hat\xi_n\theta_n}\bar p_s(0)$, *not* by $g_{st}(\theta)\bar p_s(0)$. Watch which frame the point is written in.

---

## 8. The Adjoint transformation and order of exponentials (L11) ★★★

### 8.1 Definition

To change the coordinates a twist is written in, use
$$\boxed{\mathrm{Ad}_g = \begin{bmatrix}R & \hat pR\\0&R\end{bmatrix}\in\mathbb R^{6\times6},\qquad \mathrm{Ad}_g^{-1} = \mathrm{Ad}_{g^{-1}} = \begin{bmatrix}R^T & -R^T\hat p\\0&R^T\end{bmatrix}}$$
for $g = (p, R)$.

**Lemma 2.13:** $g\hat\xi g^{-1}$ is a twist, with coordinates $\mathrm{Ad}_g\xi$. Combined with 5.4:
$$e^{(\mathrm{Ad}_g\xi)^\wedge\theta} = g\,e^{\hat\xi\theta}g^{-1}$$

**What it does geometrically:** $\mathrm{Ad}_g\xi$ is the *same physical screw axis*, moved rigidly by g.

```
     original axis (ω, through q)          moved axis (Rω, through Rq + p)
            ↑ ω                                    ↗ ω' = Rω
            │                   ── g ──▶          ╱
       q ●──┘                               q' ●╱      q' = g(q) = Rq + p
```

### 8.2 The by‑hand shortcut (avoid 6×6 matrices) ★★★

For a **revolute** twist $\xi = (-\omega\times q, \omega)$:
$$\boxed{\mathrm{Ad}_g\xi = \begin{bmatrix}-\omega'\times q'\\\omega'\end{bmatrix},\quad \omega' = R\omega,\quad q' = Rq + p}$$
For a **prismatic** twist $(v, 0)$: $\mathrm{Ad}_g\xi = (Rv, 0)$.

So **move the axis and rebuild the twist**. When the motion is simple (a planar arm, one rotated joint), just read $q'$ and $\omega'$ from the geometry with trigonometry, exactly as the official review solution does:

> **Review‑sheet example (SCARA joint 3).** After joints 1 and 2 move:
> - $\omega_3' = z$
> - $q_3' = (-l_1s_1 - l_2s_{12},\ l_1c_1 + l_2c_{12},\ 0)$
>
> So $\xi_3' = \big(l_1c_1 + l_2c_{12},\ l_1s_1+l_2s_{12},\ 0,\ 0,0,1\big)$.
> Formally this is $\xi_3' = \mathrm{Ad}_{e^{\hat\xi_1\theta_1}e^{\hat\xi_2\theta_2}}\xi_3$.

### 8.3 Why the order of the exponentials works

**Move joint 1 first.** Axis 2 now sits at $\xi_2' = \mathrm{Ad}_{e^{\hat\xi_1\theta_1}}\xi_2$. Then
$$g_{st} = e^{\hat\xi_2'\theta_2}e^{\hat\xi_1\theta_1}g_{st}(0) = \big(e^{\hat\xi_1\theta_1}e^{\hat\xi_2\theta_2}e^{-\hat\xi_1\theta_1}\big)e^{\hat\xi_1\theta_1}g_{st}(0) = e^{\hat\xi_1\theta_1}e^{\hat\xi_2\theta_2}g_{st}(0)$$
**Same answer.** The physical order doesn't matter, but if you move joint 1 first you have to use the *moved* twist.

**General rule:** the current spatial twist of joint i is
$$\boxed{\xi_i' = \mathrm{Ad}_{e^{\hat\xi_1\theta_1}\cdots e^{\hat\xi_{i-1}\theta_{i-1}}}\,\xi_i}$$
This depends only on the joints **before** i. Remember this; it is exactly column i of the spatial Jacobian (Section 13).

### 8.4 Expressing a twist in another frame (e.g., the tool frame) ★★★

If $\xi^s$ is written in S, then the same twist in frame B is $\xi^b = \mathrm{Ad}_{g_{sb}}^{-1}\xi^s$.

**The much easier way to do it:** a joint axis is **rigidly fixed to every link after it**. So the axis of joint i, written in the coordinates of any link frame $j\ge i$ (including the tool when j = n), is **constant**. You can read it straight off the reference configuration:
$$\text{(joint } n\text{'s twist in tool coords, any } \theta) = \mathrm{Ad}^{-1}_{g_{st}(\theta)}\xi_n' = \mathrm{Ad}^{-1}_{g_{st}(0)}\xi_n$$

**Algebraic proof:** $g_{st}(\theta) = (e^{1}\cdots e^{n-1})\,e^{n}\,g_{st}(0)$, and $e^{\hat\xi_n\theta_n}$ commutes with its own twist, so $\mathrm{Ad}_{e^{-\hat\xi_n\theta_n}}\xi_n = \xi_n$.

**In practice:** write ω and q of the last joint *as seen from the tool frame* in the picture, then form $(-\omega\times q, \omega)$. With frames parallel and the tool at distance $\ell$ along +y from the last axis point, for example, $q^T = (0, -\ell, 0)$.

---

## 9. Denavit–Hartenberg as relative twists (L9, L11) ★★★ (Q6b–d)

**This course's D‑H is not the textbook "a, α, d, θ" table.** It is POE written link to link:
$$g_{st}(\theta) = g_{sl_1}(\theta_1)\,g_{l_1l_2}(\theta_2)\cdots g_{l_{n-1}l_n}(\theta_n)\,g_{l_nt},\qquad \boxed{g_{l_{i-1}l_i}(\theta_i) = e^{\hat\xi_{i-1,i}\theta_i}\,g_{l_{i-1}l_i}(0)}$$

Here **$\xi_{i-1,i}$ is joint i's twist written in the coordinates of link frame $L_{i-1}$** (not S).
- Link frames sit at each link's proximal joint ($L_i$ on joint i's axis).
- $L_0 = S$, so $\xi_{0,1} = \xi_1$.

**Why it's constant:** axis i is fixed to link i−1, so its coordinates in $L_{i-1}$ never change, whatever the other joints do. You can define it at any configuration and get the same numbers.

**The relation to POE (L11):**
$$\boxed{\xi_i = \mathrm{Ad}_{g_{sl_{i-1}}(0)}\,\xi_{i-1,i}\qquad\Longleftrightarrow\qquad \xi_{i-1,i} = \mathrm{Ad}^{-1}_{g_{sl_{i-1}}(0)}\,\xi_i}$$

**By‑hand rule when the link frames are parallel to S** (the usual exam assumption):
- revolute: same ω, but measure q **from $L_{i-1}$'s origin**: $\xi_{i-1,i} = \big(-\omega_i\times(q_i - o_{i-1}),\ \omega_i\big)$
- prismatic: **identical** to the POE twist, $(v_i, 0)$ (moving a frame doesn't change a direction)
- $g_{l_nt}$ = the constant offset from the last link frame to the tool. Note $g_{l_nt}\neq g_{st}(0)$.

> **Worked example (continuing the 4‑DOF robot from 7.3).**
> Link frames, all parallel to S:
> - $L_1 = L_2 = L_3$ at $(0,0,l_0)$, where axes 1 and 2 meet. The prismatic joint's frame can sit anywhere on its line.
> - $L_4$ at $(0,l_1,l_0)$.
>
> | i | POE $\xi_i$ | D‑H $\xi_{i-1,i}$ (from $L_{i-1}$) | why |
> |---|---|---|---|
> | 1 | $(0,0,0,\ 0,0,1)$ | $(0,0,0,\ 0,0,1)$ | $L_0 = S$ |
> | 2 | $(0,l_0,0,\ 1,0,0)$ | $(0,0,0,\ 1,0,0)$ | axis 2 passes through $L_1$'s origin |
> | 3 | $(0,1,0,\ 0,0,0)$ | $(0,1,0,\ 0,0,0)$ | prismatic: unchanged |
> | 4 | $(0,l_0,-l_1,\ 1,0,0)$ | $(0,0,-l_1,\ 1,0,0)$ | $q - o_3 = (0,l_1,0)$, and $(0,l_1,0)\times x = (0,0,-l_1)$ |
>
> And $g_{l_4t} = [I \mid (0,l_2,0)]$.
>
> **Check row 4 with the Adjoint:** $\mathrm{Ad}_{g_{sl_3}(0)}$ with $p = (0,0,l_0)$ gives $v' = v + p\times\omega = (0,0,-l_1) + (0,0,l_0)\times(1,0,0) = (0,0,-l_1) + (0,l_0,0) = (0,l_0,-l_1)$ ✓.

**Relative transform, the 4×4 matrix you're asked to write out** (e.g., $g_{l_3l_4}(\theta_4)$):
$$g_{l_3l_4}(\theta_4) = \begin{bmatrix}R_x(\theta_4) & (I - R_x(\theta_4))q\\0&1\end{bmatrix}\begin{bmatrix}I & (0,l_1,0)^T\\0&1\end{bmatrix} = \begin{bmatrix}R_x(\theta_4) & (0,l_1,0)^T\\0&1\end{bmatrix},\quad q = (0,l_1,0)$$

> **Shortcut ★★:** when $L_i$ sits **on** joint i's axis, the result is always
> $$g_{l_{i-1}l_i}(\theta_i) = \begin{bmatrix}e^{\hat\omega_i\theta_i} & p_{l_{i-1}l_i}(0)\\0&1\end{bmatrix}$$
> The frame just spins in place. (Algebra: $(I-R)p + Rp = p$.)
> For a **prismatic** joint, it's $[R_{l_{i-1}l_i}(0) \mid p(0) + v\theta_i]$.

---

## 10. Inverse kinematics and the Paden–Kahan subproblems (L12–L13, HW3) ★★★

### 10.1 The problem

Given $g_d$, solve $g_{st}(\theta) = g_d$.
- Move the constant to the right: $e^{\hat\xi_1\theta_1}\cdots e^{\hat\xi_n\theta_n} = g_dg_{st}(0)^{-1} =: g_1$, which is known.
- Each subproblem can have 0, 1, 2, or ∞ solutions.
- **The maximum number of solutions belongs to the robot and $g_d$**, not to the solution method you used.

### 10.2 The two tricks that power everything ★★★

1. **A point on a revolute axis doesn't move:** $e^{\hat\xi\theta}\bar p = \bar p$ for all θ, if p is on the zero‑pitch axis of ξ. Apply both sides of the equation to such points to **kill exponentials**.
2. **Rigid motions preserve distances:** subtract a point fixed by the first joints, then take norms to **kill exponentials on the left**.

### 10.3 Warm‑up: the planar RR arm (L12)

$x = l_1c_1 + l_2c_{12}$, $y = l_1s_1 + l_2s_{12}$, $r = \sqrt{x^2+y^2}$, $\phi = \operatorname{atan2}(y,x)$.

$$\cos\theta_2 = \frac{r^2 - l_1^2 - l_2^2}{2l_1l_2}\ \Rightarrow\ \theta_2 = \pm\cos^{-1}(\cdot)\qquad \theta_1 = \phi - \operatorname{atan2}(l_2s_2,\ l_1 + l_2c_2)$$

That's 2 solutions (elbow up/down), 1 at full extension or full fold, and 0 if $r > l_1+l_2$ or $r < |l_1-l_2|$.
- Lecture form: $\theta_2 = \pi\pm\alpha$ with $\alpha = \cos^{-1}\frac{l_1^2+l_2^2-r^2}{2l_1l_2}$, and $\theta_1 = \phi\pm\beta$ with $\beta = \cos^{-1}\frac{r^2+l_1^2-l_2^2}{2l_1r}$, signs paired.

**Redundant robots:** a planar RRRR arm has 4 joints but SE(2) has only 3 DOF ($x, y, \phi$). That means **infinitely many** solutions (a self‑motion), not "at most 4." Similarly, the 7‑DOF Canadarm2 is redundant in SE(3).

### 10.4 Subproblem 1: rotation about one axis ★★★

**Solve** $e^{\hat\xi\theta}\bar p = \bar q$, with ξ zero‑pitch, unit ω, and r on the axis.

```
   project onto the plane ⟂ ω:          side view:
                                             ω ↑
        v'                                    │   ● q
        ●                                     │  ╱      p and q must be at the
         ╲  θ                                 │ ╱       SAME height along ω and
          ╲ ↺                                r●─────● p  the SAME radius from the axis
           ●────────● u'
           r(axis)
```

**Steps:**
1. $u = p - r$, $v = q - r$
2. $u' = u - \omega\omega^Tu$, $v' = v - \omega\omega^Tv$ (projections onto the plane ⟂ ω)
3. **Solvable iff** $\omega^Tu = \omega^Tv$ (same height) **and** $\|u'\| = \|v'\|$ (same radius)
4. $$\boxed{\theta = \operatorname{atan2}\big(\omega^T(u'\times v'),\ u'^Tv'\big)}$$
5. **One solution**, or **infinitely many** if $u' = 0$ (p on the axis).

> **T/F ★★:** "A solution to $e^{\hat\xi\theta}p = q$ always exists as long as $p\neq q$." **False.** You need equal height and equal radius. Most pairs of points can't be reached by one rotation.

**By hand with a coordinate axis:** for $\omega = z$ and $r = 0$, θ is just $\operatorname{atan2}(q_y, q_x) - \operatorname{atan2}(p_y, p_x)$.

### 10.5 Subproblem 2: two intersecting axes ★★★

**Solve** $e^{\hat\xi_1\theta_1}e^{\hat\xi_2\theta_2}\bar p = \bar q$. The axes meet at r.

**Idea:** find the intermediate point $c = e^{\hat\xi_2\theta_2}p = e^{-\hat\xi_1\theta_1}q$. It must lie on **both**:
- the circle p sweeps around axis 2
- the circle q sweeps around axis 1

Two circles on a sphere of radius $\|p - r\|$ intersect in 0, 1, or 2 points.

```
              ω1
              ↑
              │    ● q
              │   ╱ (circle of q about ω1)
          ● c ┤  ╱
             ╲│ ╱
       r ●────┼──────▶ ω2
               ╲
                ● p  (circle of p about ω2)
```

**Formula** ($u = p-r$, $v = q - r$, $z = c - r = \alpha\omega_1 + \beta\omega_2 + \gamma(\omega_1\times\omega_2)$):
$$\alpha = \frac{(\omega_1^T\omega_2)\,\omega_2^Tu - \omega_1^Tv}{(\omega_1^T\omega_2)^2 - 1},\qquad \beta = \frac{(\omega_1^T\omega_2)\,\omega_1^Tv - \omega_2^Tu}{(\omega_1^T\omega_2)^2-1}$$
$$\gamma^2 = \frac{\|u\|^2 - \alpha^2 - \beta^2 - 2\alpha\beta\,\omega_1^T\omega_2}{\|\omega_1\times\omega_2\|^2}$$

- $\gamma^2 > 0$ gives two solutions, $\gamma^2 = 0$ one, and $\gamma^2 < 0$ none.
- **Then do SP1 twice:** rotate p → c about $\omega_2$ to get $\theta_2$, and rotate c → q about $\omega_1$ to get $\theta_1$.

**Perpendicular axes** ($\omega_1^T\omega_2 = 0$), the common exam case:
$$\alpha = \omega_1^Tv,\qquad \beta = \omega_2^Tu,\qquad \gamma^2 = \|u\|^2 - \alpha^2 - \beta^2$$

**Intuition:** $\omega_1$ preserves q's height along $\omega_1$, and $\omega_2$ preserves p's height along $\omega_2$. So c inherits both heights.

> **Worked example (not HW3's numbers).** $\omega_1 = z$, $\omega_2 = x$, both through $r = 0$. Take $p = (0,1,0)$ and $q = (-\tfrac{1}{\sqrt2}, 0, \tfrac{1}{\sqrt2})$.
>
> 1. $\alpha = z^Tq = \tfrac{1}{\sqrt2}$, $\beta = x^Tp = 0$, $\gamma^2 = 1 - \tfrac12 - 0 = \tfrac12$, so $\gamma = \pm\tfrac{1}{\sqrt2}$.
> 2. $\omega_1\times\omega_2 = z\times x = y$, so $c = (0, \pm\tfrac{1}{\sqrt2}, \tfrac{1}{\sqrt2})$.
> 3. **θ₂** (SP1 about x, p → c): $u' = (0,1,0)$ and $v' = c$. Then $x^T(u'\times v') = \tfrac{1}{\sqrt2}$ and $u'^Tv' = \pm\tfrac{1}{\sqrt2}$, so $\theta_2 = 45°$ or $135°$.
> 4. **θ₁** (SP1 about z, c → q): $u' = (0,\pm\tfrac{1}{\sqrt2},0)$ and $v' = (-\tfrac1{\sqrt2},0,0)$. Then $z^T(u'\times v') = \pm\tfrac12$ and $u'^Tv' = 0$, so $\theta_1 = \pm90°$.
>
> **Two solutions:** $(90°, 45°)$ and $(-90°, 135°)$. Check the first: $R_x(45°)p = (0, .707, .707)$, then $R_z(90°)$ gives $(-.707, 0, .707) = q$ ✓.

### 10.6 Subproblem 3: rotate to a given distance ★★★

**Solve** $\|\bar q - e^{\hat\xi\theta}\bar p\| = \delta$.

The circle of p around the axis meets the sphere of radius δ around q. That gives **0, 1, or 2** solutions.

**Steps:**
1. $u, v, u', v'$ as in SP1.
2. $\delta'^2 = \delta^2 - |\omega^T(p - q)|^2$ (remove the height difference)
3. $\theta_0 = \operatorname{atan2}(\omega^T(u'\times v'),\ u'^Tv')$
4. $$\boxed{\theta = \theta_0\pm\cos^{-1}\!\Big(\frac{\|u'\|^2 + \|v'\|^2 - \delta'^2}{2\|u'\|\|v'\|}\Big)}$$

The last step is the law of cosines in the projected plane.

### 10.7 Full IK strategy: elbow manipulator (L12–13) ★★★

$e^1e^2e^3e^4e^5e^6 = g_1 := g_dg_{st}(0)^{-1}$, with wrist center $p_w$ (axes 4, 5, 6 meet there) and $p_b$ (axes 1, 2 meet there).

| Step | Trick | Subproblem | Number of solutions |
|---|---|---|---|
| **θ₃** | Apply to $p_w$: the wrist exponentials vanish, giving $e^1e^2e^3p_w = g_1p_w$. Subtract $p_b$ and take norms: $\|e^3p_w - p_b\| = \|g_1p_w - p_b\| = \delta$ | SP3 | 2 |
| **θ₁, θ₂** | With θ₃ known, $e^1e^2(e^3p_w) = g_1p_w$ | SP2 | 2 |
| **θ₄, θ₅** | $e^4e^5e^6 = e^{-3}e^{-2}e^{-1}g_1 =: g_2$. Apply to a point on axis 6 but not on axes 4 or 5 (for example the tool origin) | SP2 | 2 |
| **θ₆** | $e^6 = g_3$ (known). Apply to a point **not** on axis 6 | SP1 | 1 |

**Total: up to 2·2·2 = 8 solutions.**

**Planning checklist for "show how you'd solve the IK":**
1. Look for **intersecting axes** (they enable SP2) and **wrist centers** (points fixed by the last 3 joints).
2. Look for **prismatic joints**. They often control one coordinate directly. On SCARA, θ₄ = z − l₀ by inspection.
3. Find points on axes that kill exponentials. Use distance (SP3) to isolate a single joint.
4. Peel off solved joints: move their inverses to the right‑hand side.
5. **Count solutions** by multiplying the branch counts.

**SCARA (review sheet):**
- θ₄ from the z‑coordinate.
- θ₂ from SP3 (p on axis 3, q on axis 1).
- θ₁ from SP1.
- θ₃ from SP1.
- Rotation: $\theta_1+\theta_2+\theta_3 = \phi$.

**Facts to remember:**
- A **prismatic + revolute pair on the same axis commutes** (rotating about and sliding along one line don't interfere). Swapping their order gives the same IK solutions.
- Given all the other joints, a single revolute angle in $[0, 2\pi)$ is **unique**, unless the point you'd use lies on its axis.
- **At the boundary of the workspace** (fully stretched), the elbow up/down pair merges into one. Other branches (base flipped by π with the shoulder "over the top," wrist flips) can still exist, so "only a single solution" is generally false for 6‑DOF arms.

---

## 11. Rigid body velocity (L13–L14, HW3) ★★★ (Q7)

### 11.1 Rotational velocity

From $RR^T = I$, differentiating gives $\dot RR^T + R\dot R^T = 0$. So **$\dot RR^T$ is skew**, and likewise $R^T\dot R$ (Lemma 2.12).

**Spatial angular velocity** (axis seen from A) and **body angular velocity** (axis seen from B):
$$\hat\omega^s_{ab} = \dot R_{ab}R_{ab}^T,\qquad \hat\omega^b_{ab} = R_{ab}^T\dot R_{ab},\qquad \omega^s = R_{ab}\,\omega^b$$

> **T/F ★★★:** "If $R(t)\in SO(3)$, then $\dot R(t)\in so(3)$." **False.** $\dot R = \hat\omega^sR$, which is not skew in general. It's $\dot RR^T$ and $R^T\dot R$ that are skew. (The only exception is at $R = I$.)

**Velocity of a point on the body:** $v_{q_a} = \omega^s\times q_a$, or in body coordinates $v_{q_b} = \omega^b\times q_b$.

### 11.2 General (SE(3)) velocity ★★★

$$\hat V^s_{ab} = \dot g_{ab}g_{ab}^{-1},\qquad V^s_{ab} = \begin{bmatrix}v^s\\\omega^s\end{bmatrix} = \begin{bmatrix}-\dot RR^Tp + \dot p\\ (\dot RR^T)^\vee\end{bmatrix}$$
$$\hat V^b_{ab} = g_{ab}^{-1}\dot g_{ab},\qquad V^b_{ab} = \begin{bmatrix}v^b\\\omega^b\end{bmatrix} = \begin{bmatrix}R^T\dot p\\ (R^T\dot R)^\vee\end{bmatrix}$$
$$\boxed{V^s = \mathrm{Ad}_{g_{ab}}V^b}\qquad V^b_{ab} = -V^s_{ba}\ \ (\text{Lemma 2.16})$$

**What each piece *means*. This is heavily tested, so read it carefully:**

| Quantity | Physical meaning |
|---|---|
| $\omega^s$ | angular velocity of B, in A's coordinates |
| $v^s$ | ⚠ **NOT** the velocity of B's origin. It's the velocity of the (imaginary) body point that is **currently passing through A's origin**. |
| $\omega^b$ | angular velocity of B, in B's coordinates |
| $v^b$ | the velocity of **B's origin**, written in B's coordinates. **Not** zero in general. |
| $V^b$ | independent of which inertial frame you use ($V^b_{ab} = V^b_{cb}$) |

```
                 body B (rotating about the axis through c)
           ┌─────────────┐
           │      ● B origin: its velocity = ω×p + v^s   (written in A)
           │   ↻ c       │
           └─────────────┘
       ● A origin: the imaginary body point here moves with velocity v^s
```

**Point velocities** (point × twist gives a vector, so the bottom entry is 0):
$$\bar v_{q_a} = \hat V^s\bar q_a = \begin{bmatrix}\omega^s\times q_a + v^s\\0\end{bmatrix},\qquad \bar v_{q_b} = \hat V^b\bar q_b = \begin{bmatrix}\omega^b\times q_b + v^b\\0\end{bmatrix}$$

### 11.3 The fastest way to get $V^s$ for a known motion ★★★

If the body **rotates about a fixed axis** (ω, through q) at rate $\dot\theta$, then $g(t) = e^{\hat\xi\theta(t)}g(0)$, and
$$\boxed{V^s = \xi\,\dot\theta = \begin{bmatrix}-\omega\times q\\\omega\end{bmatrix}\dot\theta}$$
**The spatial velocity is just the joint twist times the joint rate.** Add translation along the axis to get a screw: $v^s = -\omega\times q\,\dot\theta + v_{\text{trans}}$.

**Integrating a constant velocity:**
- constant $V^s$: $g(t) = e^{\hat V^st}g(0)$
- constant $V^b$: $g(t) = g(0)e^{\hat V^bt}$

> **Worked example (HW3‑style, different numbers).** Body B starts with its origin at $(2,0,0)$ in A, axes parallel. It rotates about its own z‑axis at $\pi$ rad/s.
> - **Spatial velocity:** $\omega = (0,0,\pi)$, axis through $q = (2,0,0)$. $v^s = -\omega\times q = -\pi\,(z\times(2,0,0)) = -\pi(0,2,0)$, so $V^s = (0,-2\pi,0,\ 0,0,\pi)$.
> - **Point** $q_b = (1,0,0)$, so $q_a = (3,0,0)$. $v_{q_a} = \omega\times q_a + v^s = (0,3\pi,0) + (0,-2\pi,0) = (0,\pi,0)$.
>   - **Sanity check:** it's 1 unit from the spin axis and spinning at π, so its speed is π, tangent in +y ✓.
> - **Position after 0.5 s:** the body has rotated π/2 about the line x = 2, so $q_a = (2,0,0) + R_z(\tfrac\pi2)(1,0,0) = (2,1,0)$.
> - **Add translation at 3 units/s along +z:** $V^s = (0,-2\pi,3,\ 0,0,\pi)$, with pitch $h = \omega^Tv/\|\omega\|^2 = 3\pi/\pi^2 = 3/\pi$. After 0.5 s, $q_a = (2,1,1.5)$ (verified with $e^{\hat V^s\cdot0.5}$).

---

## 12. Wrenches (L14) ★★

**A wrench is a force and a moment at a point**, written in some frame's coordinates:
$$F = \begin{bmatrix}f\\\tau\end{bmatrix}\in\mathbb R^6$$
For example, $F_b$ is applied at B's origin, in B's coordinates.

**Power pairs twists with wrenches:**
$$\delta W = V^b\cdot F_b = v^b\cdot f_b + \omega^b\cdot\tau_b$$
Two wrenches are **equivalent** if they do the same work for every motion.

**Changing frames uses the transpose of the Adjoint**, not a 4×4 g:
$$\boxed{F_c = \mathrm{Ad}^T_{g_{bc}}F_b,\qquad \mathrm{Ad}_g^T = \begin{bmatrix}R^T & 0\\ -R^T\hat p & R^T\end{bmatrix}\ \Rightarrow\ f_c = R^Tf_b,\quad \tau_c = R^T(\tau_b - p\times f_b)}$$
Here p is the position of C's origin in B's coordinates.

**Where this comes from:** power is invariant. $V_c = \mathrm{Ad}^{-1}_{g_{bc}}V_b$, so $V_c^TF_c = V_b^TF_b$ forces $F_c = \mathrm{Ad}^T_{g_{bc}}F_b$.

**Physical meaning:** the force is the same vector (just rotated into new coordinates). The **moment changes** with the reference point, by $-p\times f$. So the **net wrench depends on where you measure it**: moving the reference point P changes the moment part.

> **T/F ★★:** "To move a wrench from frame A to frame B, multiply by a homogeneous transform in SE(3)." **False.** You need a 6×6 $\mathrm{Ad}^T$. (A 4×4 g acts on points, not wrenches.)

---

## 13. The manipulator Jacobian (MLS §3.4, *probably upcoming in lecture*; tested on both exams) ★★★ (Q7)

### 13.1 The spatial Jacobian

$$V^s_{st} = J^s_{st}(\theta)\,\dot\theta,\qquad \boxed{J^s_{st}(\theta) = \big[\ \xi_1\ \ \xi_2'\ \ \cdots\ \ \xi_n'\ \big],\quad \xi_i' = \mathrm{Ad}_{e^{\hat\xi_1\theta_1}\cdots e^{\hat\xi_{i-1}\theta_{i-1}}}\,\xi_i}$$

**Column i is joint i's twist, at the current configuration, in S coordinates.** You compute it with the shortcut from 8.2: move the axis point and direction using the earlier joints, then rebuild $(-\omega'\times q', \omega')$.
- Column 1 is always the reference twist $\xi_1$.
- **$J^s$ never depends on $\theta_n$.** Column i depends only on $\theta_1,\dots,\theta_{i-1}$.
- **At $\theta = 0$:** $J^s(0) = [\xi_1\ \cdots\ \xi_n]$, just the reference twists. That's why "assume θ = 0" makes exam problems easy.

### 13.2 The body Jacobian

$$V^b_{st} = J^b_{st}\dot\theta,\qquad J^b_{st} = \big[\xi_1^\dagger\ \cdots\ \xi_n^\dagger\big],\quad \xi_i^\dagger = \mathrm{Ad}^{-1}_{e^{\hat\xi_i\theta_i}\cdots e^{\hat\xi_n\theta_n}g_{st}(0)}\,\xi_i$$

**Column i is joint i's twist in tool coordinates.**
- It depends only on $\theta_{i+1},\dots,\theta_n$. The $e^{\hat\xi_i\theta_i}$ in the subscript drops out, because a twist is invariant under its own exponential.
- **So moving joint i does not change column i of $J^b$.**
- **$J^b$ never depends on $\theta_1$.** It *can* depend on $\theta_n$.
- **The last column is constant:** $\xi_n^\dagger = \mathrm{Ad}^{-1}_{g_{st}(0)}\xi_n$ (see 8.4).
- **Conversion:** $J^s = \mathrm{Ad}_{g_{st}}J^b$.

### 13.3 Worked example (3‑DOF, my own robot) ★★★

```
          z ↑   T (0,0,l1+l2)          J1: revolute about +z through 0
            │                          J2: revolute about +x through 0
            ●  J3 (0,0,l1), axis +x    J3: revolute about +x through (0,0,l1)
            │                          Arm points straight up at θ = 0;
            ●  J1/J2 at S              frames parallel.
            └──────▶ y
```

- **Reference twists:**
  - $\xi_1 = (0,0,0,\ 0,0,1)$
  - $\xi_2 = (0,0,0,\ 1,0,0)$
  - $\xi_3$: $q\times\omega = (0,0,l_1)\times(1,0,0) = (0, l_1, 0)$, so $\xi_3 = (0,l_1,0,\ 1,0,0)$
  - $g_{st}(0) = [I\mid(0,0,l_1+l_2)]$
- **Column 2:** joint 1 rotates axis 2 by $R_z(\theta_1)$. So $\omega_2' = (c_1, s_1, 0)$, with $q_2' = 0$, which gives $\xi_2' = (0,0,0,\ c_1,s_1,0)$.
- **Column 3:**
  - $\omega_3' = R_z(\theta_1)R_x(\theta_2)\,x = (c_1,s_1,0)$
  - $q_3' = R_z(\theta_1)R_x(\theta_2)(0,0,l_1) = R_z(\theta_1)(0,-l_1s_2,l_1c_2) = (l_1s_1s_2,\ -l_1c_1s_2,\ l_1c_2)$
  - $v_3' = q_3'\times\omega_3' = (-l_1s_1c_2,\ l_1c_1c_2,\ l_1s_2)$

$$J^s_{st}(\theta) = \begin{bmatrix}0&0&-l_1s_1c_2\\0&0&l_1c_1c_2\\0&0&l_1s_2\\0&c_1&c_1\\0&s_1&s_1\\1&0&0\end{bmatrix}\qquad(\text{checked by differentiating } g_{st}(\theta)\text{ symbolically})$$

**Velocity of a tool point at θ = 0** ($\bar q^s = (q_x,q_y,q_z,1)$):
- $V^s = J^s(0)\dot\theta = (0,\ l_1\dot\theta_3,\ 0,\ \dot\theta_2+\dot\theta_3,\ 0,\ \dot\theta_1)$
- $v_q = \omega\times q + v = \big(-\dot\theta_1q_y,\ \ \dot\theta_1q_x - (\dot\theta_2+\dot\theta_3)q_z + l_1\dot\theta_3,\ \ (\dot\theta_2+\dot\theta_3)q_y\big)$

**Body Jacobian at θ = 0:** $\xi_i^\dagger = \mathrm{Ad}^{-1}_{g_{st}(0)}\xi_i$. That's just each axis re‑measured from the tool origin $(0,0,l_1+l_2)$:
- $\xi_1^\dagger = (0,0,0,\ 0,0,1)$
- $\xi_2^\dagger$: axis through $(0,0,-(l_1+l_2))$ in tool coordinates, so $\xi_2^\dagger = (0,-(l_1+l_2),0,\ 1,0,0)$
- $\xi_3^\dagger = (0,-l_2,0,\ 1,0,0)$

**Joint torques for a body wrench $F_b = (f_x,f_y,f_z,\tau_x,\tau_y,\tau_z)$ at the tool, θ = 0:**
$$\tau = (J^b_{st})^TF_b = \big(\tau_z,\ \ \tau_x - (l_1+l_2)f_y,\ \ \tau_x - l_2f_y\big)$$
**Sanity check:** a +y force at height L about the x‑axis gives moment $(0,0,L)\times(0,f_y,0) = (-Lf_y, 0, 0)$ ✓.

### 13.4 Statics: $\tau = J^TF$ ★★★

$$\boxed{\tau = (J^b_{st})^TF_t\ \ (\text{tool wrench in tool coords})\qquad\text{or}\qquad \tau = (J^s_{st})^TF_s\ \ (\text{spatial wrench})}$$

**Where it comes from:** power balance, $\tau^T\dot\theta = F^TV = F^TJ\dot\theta$ for every $\dot\theta$.
- Note the **transpose**: J is 6×n and τ is n×1.
- **A wrench on link k** (not the tool) uses that link's Jacobian, $J^b_{sl_k}$. Its columns for joints after k are **zero**, because those joints can't push on link k.
- **Sign:** τ = JᵀF gives the torques for the robot to *apply* F. To *hold against* an external wrench, flip the sign. Read the wording.

### 13.5 Singularities ★★

**A configuration is singular when J drops rank.** The robot loses an instantaneous direction of motion. Classic ones (MLS):
1. **Two revolute axes collinear.** Example: a wrist with joint 5 at 0 lines up axes 4 and 6. That's the review sheet's answer.
2. Three parallel, coplanar revolute axes (e.g., the elbow fully stretched)
3. Four revolute axes intersecting at one point
4. Four coplanar revolute axes
5. Six revolute axes intersecting one common line
6. A prismatic axis ⟂ to the plane of two parallel coplanar revolute axes

---

## 14. True/False drill: the concepts behind every exam statement ★★★

Cover the answer column and test yourself. Each line is a concept the exams turn into a T/F question. Answers are for the statement as written here.

| # | Statement | T/F | One‑line reason |
|---|---|---|---|
| 1 | $R_{ca} = R_{bc}^TR_{ba}$ | T | $R_{cb}R_{ba}$, and $R_{cb} = R_{bc}^T$ |
| 2 | so(3) is an abelian group under matrix addition | T | it's a vector space |
| 3 | so(3) is a group under matrix multiplication | F | not closed (the product of skews isn't skew) |
| 4 | SO(2) is abelian; SO(3) is not | T | one axis in the plane vs. many axes in 3D |
| 5 | $\det\hat\omega = 1$ for unit ω | F | 3×3 skew means det = 0 ($\hat\omega\omega = 0$) |
| 6 | Rodrigues holds for every $\hat\omega\in so(3)$ and θ | F | needs $\|\omega\|=1$ |
| 7 | exp: so(3)→SO(3) is surjective | T | Prop 2.5 |
| 8 | exp: so(3)→SO(3) is injective | F | $(\omega,\theta)$ and $(-\omega, 2\pi-\theta)$ give the same R; at R = I every ω works |
| 9 | Some R has infinitely many exponential coordinates with θ∈[0,2π) | T | R = I |
| 10 | θ = π is a singularity of the axis‑angle map | F | ω is still recoverable from $\omega\omega^T = (R+I)/2$; only θ = 0 is singular |
| 11 | Orthogonal A with eigenvalues $\{1, e^{\pm i\phi}\}$ is in SO(3) | T | det = product of the eigenvalues = 1 |
| 12 | $\dot R\in so(3)$ | F | $\dot RR^T$ and $R^T\dot R$ are; $\dot R$ itself isn't |
| 13 | Euler angles → SO(3) is injective / global | F | singular at β = 0 (ZYZ); every 3‑parameter chart has singularities |
| 14 | ZYX (yaw‑pitch‑roll) is singular at pitch = ±π/2 | T | gimbal lock |
| 15 | The bottom‑right entry of $\hat\xi$ is 1 | F | it's 0 (the bottom row of a twist is all zeros) |
| 16 | $\hat\xi\bar p$ has bottom entry 0 (it's a vector) | T | velocity of a point |
| 17 | Every twist with ω ≠ 0 has $v = -\omega\times q$ for some q on the axis | F | needs $+h\omega$ unless pitch = 0 |
| 18 | A revolute twist whose axis passes through the origin has v = 0 | T | $q = 0$ |
| 19 | A unit twist has $\|\xi\|_2 = 1$ | F | only $\|\omega\| = 1$ (or $\|v\| = 1$) |
| 20 | The ω ≠ 0 exp formula also works for pure translation | F | translation gives $[I\ v\theta]$; the formula would give $I$ |
| 21 | $\det g = 0$ for a pure translation g | F | $\det g = \det R = 1$ |
| 22 | exp: se(3)→SE(3) is surjective | T | Prop 2.9 |
| 23 | For $\xi\theta$ with ξ a unit twist, the screw magnitude is θ (well, $\lvert\theta\rvert$) | T* | $M = \|\omega\theta\|$. *Watch for sign wording.* |
| 24 | $e^{\hat\xi_1}e^{\hat\xi_2} = e^{\hat\xi_1+\hat\xi_2}$ | F | twists don't commute |
| 25 | POE requires the joints to physically move last‑to‑first | F | that's only how the math is ordered; the final pose doesn't depend on order |
| 26 | $\xi_3' = \mathrm{Ad}_{e^{\hat\xi_1\theta_1}e^{\hat\xi_2\theta_2}}\xi_3$ is joint 3's current spatial twist | T | Section 8.3 |
| 27 | $e^{\hat\xi_i\theta_i}g = g$ for joint i's link frame g | F | exp fixes only points on its axis |
| 28 | $g_{st}(\theta)\bar p$ moves a point given in **S coordinates** at the reference configuration | F | $g_{st}$ takes **tool** coordinates; use $e^{1}\cdots e^{n}\bar p$ |
| 29 | $\|g(\bar p - \bar q)\| = \|\bar p - \bar q\|$ for any rigid g | T | distance preserved |
| 30 | $g_{l_0l_1}(\theta) = e^{\hat\xi_1\theta_1}g_{l_0l_1}(0)$ with the standard POE twist $\xi_1$ | T | $L_0 = S$, so $\xi_{0,1} = \xi_1$ |
| 31 | D‑H twist $\xi_{i-1,i}$ is the same at any configuration (in $L_{i-1}$) | T | axis i is fixed to link i−1 |
| 32 | The max number of IK solutions depends on the solution method | F | it's a property of the robot and $g_d$ |
| 33 | A planar RRRR arm has ≤ 4 IK solutions for $g_d\in SE(2)$ | F | redundant, so infinitely many |
| 34 | A point on a zero‑pitch axis satisfies $e^{\hat\xi\theta}p = p$ for all θ | T | it's on the axis |
| 35 | $e^{\hat\xi\theta}p = q$ is solvable whenever p ≠ q | F | needs equal height and equal radius |
| 36 | A fully stretched 6R arm has exactly one IK solution | F | the elbow pair merges, but base and wrist flips remain |
| 37 | Swapping a revolute and a prismatic joint on the same axis changes the IK solutions | F | they commute |
| 38 | Every $V^s$ and $V^b$ has an associated screw (pitch 0, finite, or ∞) | T | twists ↔ screws (Chasles) |
| 39 | $v^s$ is the velocity of the body frame's origin | F | it's the body point at **A's** origin |
| 40 | The velocity of B's origin in body coordinates ($v^b$) can be non‑zero | T | $\hat V^b\bar 0 = (v^b, 0)$ |
| 41 | A wrench changes frame through a 4×4 homogeneous transform | F | it uses $\mathrm{Ad}^T$ (6×6) |
| 42 | The net wrench is independent of the reference point | F | the moment shifts by $p\times f$ |
| 43 | Spatial Jacobian $J^s$ doesn't depend on $\theta_n$ | T | column i uses only $\theta_{1..i-1}$ |
| 44 | Body Jacobian $J^b$ doesn't depend on $\theta_n$ | F | it doesn't depend on $\theta_1$ instead |
| 45 | Rotating joint i changes column i of $J^b$ | F | column i depends only on $\theta_{i+1..n}$ |
| 46 | $\tau = J^b_{st}F$ resists a wrench on link 3 | F | needs the **transpose**, and **link 3's** Jacobian |
| 47 | At θ = 0, $J^s = [\xi_1\cdots\xi_n]$ | T | all the Adjoints are identity |

---

## 15. Exam‑day recipes and common mistakes

### Q6‑type (forward kinematics, 20 points)
1. **Mark S and T on the figure.** Note x out of the page, and any odd conventions the problem forces on you.
2. **$g_{st}(0)$:** write the tool origin from the dimension lines, and R from the tool axes (usually I).
3. **For each joint:** decide revolute or prismatic. Get ω from the arrow plus the right‑hand rule. Pick the q with the most zeros. Write $\xi = (q\times\omega, \omega)$ or $(v,0)$.
4. **Check** $v\perp\omega$, and that each prismatic twist has $\omega = 0$.
5. **D‑H:** re‑measure each q from the previous link frame's origin. Prismatic twists don't change. $\xi_{0,1} = \xi_1$.
6. **Relative transform:** $g_{l_{i-1}l_i}(\theta_i) = e^{\hat\xi_{i-1,i}\theta_i}g_{l_{i-1}l_i}(0)$, multiplied out to 4×4. If $L_i$ is on axis i, it's $[e^{\hat\omega\theta}\mid p(0)]$.
7. **"Twist in tool coordinates at θ ≠ 0":** for the last joint, it's constant. Read it from the picture as seen from T. Otherwise use $\mathrm{Ad}^{-1}_{g_{st}(\theta)}\xi_i'$.
8. **Write the reasoning in one line per joint**, e.g. "ω₃ = −x (right‑hand rule on drawn arrow), q₃ = (0, ℓ₁, ℓ₀)." No reasoning means no points.

### Q7‑type (velocity, Jacobian, wrench, 20 points)
1. **Twists first** (they're partial credit by themselves).
2. **$J^s$:** column 1 is ξ₁. For each later column, rotate and translate ω and q by the preceding joints, then $(q'\times\omega', \omega')$.
3. **"Velocity of a point when θ = 0":** $V^s = [\xi_1\cdots\xi_n]\dot\theta$, then $v_q = \omega\times q + v$.
4. **"Spatial velocity when only joint k moves at rate r":** $V^s = \xi_k'r$ (its current twist times r).
5. **Torques:** $\tau = J^TF$. Match frames: body wrench goes with $J^b$, spatial wrench with $J^s$. At θ = 0, $J^b(0)$ is the twists re‑measured from the tool origin.
6. **$g_{bc}$ between two links:** $[R(\text{joint angle}) \mid \text{offset}]$, where the offset is the link length along the direction it points at reference. Example: $(0,0,-l)$ for a link hanging down.

### Mistakes that cost points
- **Sign of ω:** the drawn arrow wins over "natural" axes.
- **Mixing up $(v,\omega)$ order.** MLS puts **v on top**.
- **Writing $-q\times\omega$** (the wrong sign). The correct form is $v = -\omega\times q = q\times\omega$.
- **Forgetting a prismatic twist is $(v, 0)$**, not $(0, v)$.
- **Using Rodrigues with a non‑unit ω.**
- **Using $g_{st}(0)$ in place of $g_{l_nt}$ in D‑H**, or vice versa.
- **$\tau = JF$** without the transpose.
- **Treating $v^s$ as the origin's velocity.**
- **Skipping the fixed‑point check:** a point on a revolute axis must map to itself under $e^{\hat\xi\theta}$. Use it to catch errors.

---

## 16. What's next

1. **Today:** read this guide. Mark anything that felt shaky with a ❓.
2. **This week, light:** 10–15 minutes a day on hand computations:
   - write the twists for the four HW2 robots from memory
   - one SP1/SP2 by hand
   - one $J^s$ column by hand
3. **From Oct 26:** full practice exams in this exact format (4–5 T/F blocks plus a 20‑point FK problem plus a 20‑point Jacobian/velocity/wrench problem), all doable by hand. Save the two real midterms for timed runs in the last few days.
