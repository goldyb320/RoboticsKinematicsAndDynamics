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

**What that means for your time.** Half the exam is one long problem on POE twists, and one on Jacobians and velocities. The other 30 points are conceptual T/F questions that hinge on knowing *exactly* what a definition says (for example, "is $\dot R$ in so(3)?", "does Rodrigues' formula need $\Vert\omega\Vert=1$?").

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

1. **Preserve length:** $\Vert g(p)-g(q)\Vert= \Vert p-q\Vert$ for all points.
2. **Preserve the cross product:** $g_{\ast}(v\times w) = g_{\ast}(v)\times g_{\ast}(w)$ for all vectors.

Rule 2 is what excludes **reflections**. A mirror preserves lengths but flips handedness, so $x\times y$ would map to $-z$ instead of $z$.

Consequence: inner products are preserved, $v^Tw = g_{\ast}(v)^Tg_{\ast}(w)$, so orthogonal vectors stay orthogonal.

**Points vs. vectors.** A vector is the difference of two points, $v = q-p$, and $g_{\ast}(v) = g(q)-g(p)$. For a general motion $g=(p,R)$, a vector only gets **rotated**: $g_{\ast}(v) = Rv$. Translation cancels out in the difference. This matters later for homogeneous coordinates.

### 1.2 Rotation matrices: what the columns mean ★★★

$$
R_{ab} = \begin{bmatrix} x_{ab} & y_{ab} & z_{ab}\end{bmatrix}
$$

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

$$
R_{ac} = R_{ab}R_{bc},\qquad R_{ba} = R_{ab}^{-1} = R_{ab}^T
$$

> **T/F pattern ★★★:** "A has orientation $R_1$ relative to B, and C has orientation $R_2$ relative to B. What is A relative to C?" Translate the words into subscripts: $R_1 = R_{BA}$ and $R_2 = R_{BC}$. Then $R_{CA} = R_{CB}R_{BA} = R_2^TR_1$. Always convert words into subscripts and chain them so the inner letters match.

### 1.3 SO(3), groups, and abelian groups ★★

$$
SO(3) = \lbrace R\in\mathbb R^{3\times3} : RR^T = R^TR = I,\ \det R = +1\rbrace
$$

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

$$
a\times b = \hat a\, b,\qquad \hat a = \begin{bmatrix}0&-a_3&a_2\cr  a_3&0&-a_1\cr  -a_2&a_1&0\end{bmatrix},\qquad \hat a^T = -\hat a
$$

**Memory trick:** zeros on the diagonal; $+a_1$ sits at position (3,2), $+a_2$ at (1,3), $+a_3$ at (2,1), and each mirror entry is the negative.

**Vee map:** $(\hat a)^\vee = a$. You read $a$ back off as $(\hat a_{32}, \hat a_{13}, \hat a_{21})$.

**Key facts:**
- $\widehat{Rw} = R\,\hat w\,R^T$. Rotating the axis is the same as conjugating the skew matrix.
- $\hat w\,w = w\times w = 0$, so **$w$ is in the null space of $\hat w$**.
- **$\hat w$ is never invertible.** Its determinant is always 0, because $\hat w w = 0$ and its eigenvalues are $0, \pm i\Vert w\Vert$.

> **Trap ★★:** "If $\Vert\omega\Vert=1$, then $\hat\omega$ is invertible with $\det\hat\omega = 1$." **False.** A 3×3 skew matrix always has $\det = 0$. (The determinant of an odd‑dimension skew matrix is $\det(\hat w) = \det(\hat w^T) = \det(-\hat w) = (-1)^3\det\hat w$, so it must be 0.)

---

## 2. Euler angles (L3–L4) ★★

### 2.1 Elementary rotations (memorize these)

$$
R_x(\phi)=\begin{bmatrix}1&0&0\cr 0&c&-s\cr 0&s&c\end{bmatrix}\quad R_y(\beta)=\begin{bmatrix}c&0&s\cr 0&1&0\cr -s&0&c\end{bmatrix}\quad R_z(\alpha)=\begin{bmatrix}c&-s&0\cr s&c&0\cr 0&0&1\end{bmatrix}
$$

**Pattern:** the axis you rotate about gets the "1" row and column. The other 2×2 block is a 2D rotation with $-s$ in the upper right, **except $R_y$, where the $-s$ is in the lower left**. That's because the cyclic order is $z\to x$ for a y‑rotation.

### 2.2 ZYZ and ZYX (yaw–pitch–roll)

**Intrinsic** (rotate about the body's *new* axes each time) means you **post‑multiply** in the order you say them:

- **ZYZ:** rotate about z by α, then the new y by β, then the newest z by γ:

$$
R_{ab} = R_z(\alpha)R_y(\beta)R_z(\gamma)=\begin{bmatrix}c_\alpha c_\beta c_\gamma - s_\alpha s_\gamma & -c_\alpha c_\beta s_\gamma - s_\alpha c_\gamma & c_\alpha s_\beta\cr  s_\alpha c_\beta c_\gamma + c_\alpha s_\gamma & -s_\alpha c_\beta s_\gamma + c_\alpha c_\gamma & s_\alpha s_\beta\cr  -s_\beta c_\gamma & s_\beta s_\gamma & c_\beta\end{bmatrix}
$$

- **ZYX (yaw φ, pitch θ, roll ψ):** $R = R_z(\phi)R_y(\theta)R_x(\psi)$. Its bottom row is $(-s_\theta,\ c_\theta s_\psi,\ c_\theta c_\psi)$ and its first column is $(c_\phi c_\theta, s_\phi c_\theta, -s_\theta)$.

**Extrinsic** (about the fixed world axes) means you **pre‑multiply**. The same three rotations done about fixed axes give the product in reverse order. The intrinsic sequence z, y′, z″ equals the extrinsic sequence z, y, z with the matrices multiplied the other way around.

### 2.3 Inverse: from R back to angles

**ZYZ**, when $\sin\beta\neq0$:

$$
\beta = \mathrm{atan2}\!\big(\sqrt{r_{31}^2+r_{32}^2},\ r_{33}\big),\quad \alpha = \mathrm{atan2}(r_{23}/s_\beta,\ r_{13}/s_\beta),\quad \gamma = \mathrm{atan2}(r_{32}/s_\beta,\ -r_{31}/s_\beta)
$$

**Where these come from:** look at the third column, $(c_\alpha s_\beta, s_\alpha s_\beta, c_\beta)$, and the third row, $(-s_\beta c_\gamma, s_\beta s_\gamma, c_\beta)$. You don't need to memorize the formulas, just how to read them off the matrix.

**ZYX:** $\theta = \mathrm{atan2}(-r_{31}, \sqrt{r_{11}^2+r_{21}^2})$, $\phi = \mathrm{atan2}(r_{21}, r_{11})$, $\psi = \mathrm{atan2}(r_{32}, r_{33})$.

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

$$
\dot q(t) = \omega\times q(t) = \hat\omega\,q(t).
$$

This is a linear time‑invariant ODE, so its solution is $q(t) = e^{\hat\omega t}q(0)$. After rotating by angle $\theta$:

$$
\boxed{R(\omega,\theta) = e^{\hat\omega\theta}}
$$

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

$so(3) = \lbrace S\in\mathbb R^{3\times3}: S^T = -S\rbrace$, a 3‑dimensional real vector space.

**Lemma 2.3:**

$$
\hat a^2 = aa^T - \Vert a\Vert^2I,\qquad \hat a^3 = -\Vert a\Vert^2\hat a
$$

With $a = \omega\theta$ and $\Vert\omega\Vert=1$, the powers cycle:

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

$$
e^{\hat\omega\theta} = I + \underbrace{\Big(\theta - \tfrac{\theta^3}{3!} + \cdots\Big)}_{\sin\theta}\hat\omega + \underbrace{\Big(\tfrac{\theta^2}{2!} - \tfrac{\theta^4}{4!} + \cdots\Big)}_{1-\cos\theta}\hat\omega^2
$$

$$
\boxed{e^{\hat\omega\theta} = I + \hat\omega\sin\theta + \hat\omega^2(1-\cos\theta)}\qquad(\Vert\omega\Vert= 1\ \text{required!})
$$

> **T/F ★★★:** "Rodrigues' formula does not hold for some $\hat\omega\in so(3)$ and $\theta$." Rodrigues' formula **only holds when $\Vert\omega\Vert=1$**. The derivation used $\hat\omega^3 = -\hat\omega$, which needs a unit vector. For a non‑unit $\omega$, first normalize: $\omega\theta = \frac{\omega}{\Vert\omega\Vert}(\Vert\omega\Vert\theta)$.

**By‑hand tip:** for a coordinate axis you never need Rodrigues. $e^{\hat z\theta} = R_z(\theta)$, $e^{\hat x\theta} = R_x(\theta)$, $e^{\hat y\theta} = R_y(\theta)$. Use Rodrigues only for tilted axes. Also useful: $\hat\omega^2 = \omega\omega^T - I$.

### 3.4 exp lands in SO(3) (Prop 2.4) ★

- **Orthogonal:** $(e^{\hat\omega\theta})^{-1} = e^{-\hat\omega\theta} = e^{\hat\omega^T\theta} = (e^{\hat\omega\theta})^T$.
- **$\det = +1$:** orthogonality gives $\det = \pm1$. At $\theta=0$, $\det I = 1$, and $\det$ is continuous in $\theta$, so it can't jump to $-1$.

### 3.5 Surjectivity: every R has an axis‑angle (Prop 2.5, Euler's Thm 2.6) ★★★

Write out Rodrigues' formula entrywise, with $v_\theta = 1-\cos\theta$:

$$
R = \begin{bmatrix}\omega_1^2v_\theta + c_\theta & \omega_1\omega_2v_\theta - \omega_3s_\theta & \omega_1\omega_3v_\theta+\omega_2s_\theta\cr  \omega_1\omega_2v_\theta+\omega_3s_\theta & \omega_2^2v_\theta+c_\theta & \omega_2\omega_3v_\theta - \omega_1s_\theta\cr  \omega_1\omega_3v_\theta - \omega_2s_\theta & \omega_2\omega_3v_\theta+\omega_1s_\theta & \omega_3^2v_\theta+c_\theta\end{bmatrix}
$$

**Step 1, the angle from the trace:**

$$
\mathrm{tr}R = 1+2\cos\theta\ \Rightarrow\ \boxed{\theta = \cos^{-1}\!\Big(\frac{\mathrm{tr}R-1}{2}\Big)}
$$

**Step 2, the axis from the skew part:**

$$
\boxed{\omega = \frac{1}{2\sin\theta}\begin{bmatrix}r_{32}-r_{23}\cr  r_{13}-r_{31}\cr  r_{21}-r_{12}\end{bmatrix}}\quad(\sin\theta\neq0)
$$

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

- $SO(2) = \lbrace R(\theta)\rbrace\cong S^1$, the circle: one parameter with period 2π.
- $\hat\omega = \begin{bmatrix}0&-\omega\cr \omega&0\end{bmatrix}$ satisfies $\hat\omega^2 = -\omega^2I$, so $e^{\hat\omega\theta} = \cos(\omega\theta)I + \sin(\omega\theta)\hat\omega/\omega = R(\omega\theta)$.
- **exp: so(2)→SO(2) is surjective but not injective** ($\omega\theta$ and $\omega\theta+2\pi$ give the same R).
- **SO(2) commutes:** $R\hat\omega R^T = \hat\omega$, $Re^{\hat\omega\theta}R^T = e^{\hat\omega\theta}$, and $\frac{d}{dt}e^{\hat\omega\theta} = \hat\omega\dot\theta e^{\hat\omega\theta} = e^{\hat\omega\theta}\hat\omega\dot\theta$.
- **In 3D this fails.** $R\hat\omega R^T = \widehat{R\omega}$, which is a *different* axis.

---

## 4. SE(3) and homogeneous coordinates (L5–L6) ★★★

### 4.1 Configuration of a body

The configuration of frame B relative to A is $g_{ab} = (p_{ab}, R_{ab})$:
- $p_{ab}$ is the position of B's origin, in A's coordinates.
- $R_{ab}$ is B's axes, in A's coordinates.

$$
SE(3) = \mathbb R^3\times SO(3),\qquad q_a = p_{ab} + R_{ab}q_b
$$

**Elements of SE(3) are RBTs (Prop 2.7):**
- Lengths: $\Vert g(p_1) - g(p_2)\Vert= \Vert R(p_1-p_2)\Vert$.
- Cross products: $g_{\ast}(v) = Rv$ for vectors.

### 4.2 Homogeneous form (this is what you'll compute with)

$$
\bar g_{ab} = \begin{bmatrix}R_{ab}&p_{ab}\cr 0&1\end{bmatrix}\in\mathbb R^{4\times4},\qquad \bar q = \begin{bmatrix}q\cr 1\end{bmatrix}\ (\text{point}),\qquad \bar v = \begin{bmatrix}v\cr 0\end{bmatrix}\ (\text{vector})
$$

**The bottom entry is the type check:**

| Expression | Bottom entry | Result |
|---|---|---|
| vector ± vector | 0 | vector |
| point + vector | 1 | point |
| point − point | 0 | vector |
| point + point | 2 | **meaningless** |

A vector ignores the translation: $\bar g\bar v = (Rv, 0)$.

**Composition and inverse:**

$$
\bar g_{ac} = \bar g_{ab}\bar g_{bc} = \begin{bmatrix}R_{ab}R_{bc} & R_{ab}p_{bc}+p_{ab}\cr 0&1\end{bmatrix},\qquad \boxed{\bar g^{-1} = \begin{bmatrix}R^T & -R^Tp\cr 0&1\end{bmatrix}}
$$

**$\det\bar g = \det R = 1$, always**, even for pure translation. This is block upper triangular, so the determinant is $\det R\cdot1$.

> **T/F ★★:** "For pure translation g, $\det g = 0$." **False.** It's 1. (Don't confuse g with the twist $\hat\xi$, which *is* singular, with $\det\hat\xi = 0$.)

---

## 5. Twists and the exponential map on SE(3) (L6–L7) ★★★

### 5.1 Twist for a rotation about an *offset* axis

Let ω be a unit axis and **q any point on the axis**. A point $p$ spinning around that axis moves with velocity

$$
\dot p = \omega\times(p-q) = \hat\omega p - \hat\omega q = \hat\omega p + v,\qquad v := -\omega\times q = q\times\omega
$$

In homogeneous form this is linear: $\dot{\bar p} = \hat\xi\,\bar p$, so $\bar p(\theta) = e^{\hat\xi\theta}\bar p(0)$, where

$$
\hat\xi = \begin{bmatrix}\hat\omega & v\cr 0&0\end{bmatrix}\in se(3)\ (4\times4),\qquad \xi = \begin{bmatrix}v\cr \omega\end{bmatrix}\in\mathbb R^6\ \text{(twist coordinates; } v \text{ on top in MLS)}
$$

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

$$
\hat\xi\bar p = \begin{bmatrix}\omega\times p + v\cr 0\end{bmatrix}
$$

That bottom 0 marks it as a vector.

### 5.2 Twist types ★★★ (know these cold)

| Joint / motion | Twist coordinates $\xi = (v,\omega)$ | Notes |
|---|---|---|
| **Revolute** (pure rotation, pitch 0) | $\begin{bmatrix}-\omega\times q\cr \omega\end{bmatrix}$, $\Vert\omega\Vert= 1$, q on the axis | the axis through the origin gives $v=0$ |
| **Prismatic** (pure translation, pitch ∞) | $\begin{bmatrix}v\cr 0\end{bmatrix}$, $\Vert v\Vert= 1$ | **where** the slider sits doesn't matter, only its direction |
| **Screw/helical** (pitch h) | $\begin{bmatrix}-\omega\times q + h\omega\cr \omega\end{bmatrix}$ | translate $h\theta$ along ω while rotating θ |

**Fast cross products you'll use constantly** (right‑hand cyclic $x\to y\to z\to x$):
- $z\times x = y$, $x\times y = z$, $y\times z = x$
- reversed order flips the sign: $x\times z = -y$, and so on
- $-\omega\times q = q\times\omega$, often easier to remember

> **Worked mini‑example.** Revolute about +x through $q = (0, a, b)$:
> $v = q\times\omega = (0,a,b)\times(1,0,0) = (a\cdot0 - b\cdot0,\ b\cdot1 - 0\cdot0,\ 0\cdot0 - a\cdot1) = (0, b, -a)$.
> So $\xi = (0, b, -a, 1, 0, 0)$.
> Notice the x‑coordinate of q never matters (it's along the axis), as expected.

> **T/F ★★:** "A unit twist has 2‑norm 1." **False** in general. *Unit* means $\Vert\omega\Vert= 1$ (or $\Vert v\Vert=1$ if $\omega=0$). The full vector has $\Vert\xi\Vert= \sqrt{\Vert v\Vert^2+1}$, which is 1 only when the axis passes through the origin.

### 5.3 The exponential map formulas ★★★

**Pure translation ($\omega = 0$).** Since $\hat\xi^2 = 0$, the series stops:

$$
e^{\hat\xi\theta} = \begin{bmatrix}I & v\theta\cr 0&1\end{bmatrix}
$$

**Rotation or screw ($\omega\neq0$, $\Vert\omega\Vert=1$):**

$$
\boxed{e^{\hat\xi\theta} = \begin{bmatrix}e^{\hat\omega\theta} & (I - e^{\hat\omega\theta})(\omega\times v) + \omega\omega^Tv\,\theta\cr 0&1\end{bmatrix}}
$$

**The by‑hand shortcut for revolute joints** (the one you actually want on the exam):

$$
\boxed{\text{pure rotation about an axis through } q:\quad e^{\hat\xi\theta} = \begin{bmatrix}e^{\hat\omega\theta} & (I - e^{\hat\omega\theta})\,q\cr 0&1\end{bmatrix}}
$$

**Why the shortcut is the same thing:** with $v = -\omega\times q$, we get $\omega\times v = q - (\omega^Tq)\omega = q_\perp$, the part of q perpendicular to the axis. Also $\omega^Tv = 0$. And $(I-R)q = (I-R)q_\perp$, because R leaves the along‑axis part alone. Pick the **simplest q on the axis** and you're done in one line.

**Geometric reading:** rotate by R about the origin, then shift by $(I-R)q$ to put the pivot back at q:

$$
p \mapsto R(p - q) + q
$$

**Screw ($h\neq0$):** add $h\theta\,\omega$ to the top‑right entry:

$$
e^{\hat\xi\theta} = \begin{bmatrix}R & (I - R)q + h\theta\omega\cr 0&1\end{bmatrix}
$$

> **T/F ★★★:** "The general formula $e^{\hat\xi\theta} = \begin{bmatrix}e^{\hat\omega\theta}&(I-e^{\hat\omega\theta})(\omega\times v)+\omega\omega^Tv\theta\cr 0&1\end{bmatrix}$ works for any twist, including pure translation." **False.** Set $\omega = 0$ and the top‑right entry is $0$, but the correct answer is $v\theta$. The formula also assumed $\Vert\omega\Vert= 1$. Translation needs its own case.

> **Worked example (L10, SCARA joint 2).** $\omega = z$, $q = (0, l_1, 0)$:
>
> $$
> e^{\hat\xi_2\theta_2} = \begin{bmatrix}R_z(\theta_2) & (I - R_z(\theta_2))\begin{bmatrix}0\cr l_1\cr 0\end{bmatrix}\cr 0&1\end{bmatrix} = \begin{bmatrix}c_2&-s_2&0& l_1s_2\cr  s_2&c_2&0& l_1(1-c_2)\cr 0&0&1&0\cr 0&0&0&1\end{bmatrix}
> $$
>
> **Sanity check:** a point on the axis, $(0, l_1, 0)$, must stay put. $R_z(0,l_1,0) + (l_1s_2, l_1-l_1c_2, 0) = (-l_1s_2 + l_1s_2,\ l_1c_2 + l_1 - l_1c_2,\ 0) = (0, l_1, 0)$ ✓.
> **Always do this check.** It catches sign errors in seconds.

### 5.4 How the formula is proven (L7, the similarity trick) ★

1. Move the frame to the closest point on the axis using $g = \begin{bmatrix}I & \omega\times v\cr 0&1\end{bmatrix}$.
2. In the new coordinates, the twist becomes $\hat\xi' = g^{-1}\hat\xi g = \begin{bmatrix}\hat\omega & h\omega\cr 0&0\end{bmatrix}$, an axis through the origin. Its exponential is easy: $\begin{bmatrix}e^{\hat\omega\theta} & h\omega\theta\cr 0&1\end{bmatrix}$.
3. Transform back: $e^{\hat\xi\theta} = g\,e^{\hat\xi'\theta}g^{-1}$.

The idea worth keeping: **$e^{g\hat\xi g^{-1}\theta} = g\,e^{\hat\xi\theta}g^{-1}$.** This identity is the whole reason the Adjoint works (Section 8).

### 5.5 Surjectivity onto SE(3) (Prop 2.9) ★★

**Given any $g = (R, p)$, there is a twist and θ with $g = e^{\hat\xi\theta}$:**
- **$R = I$** (pure translation): $v = p/\Vert p\Vert$, $\theta = \Vert p\Vert$.
- **$R\neq I$:** first solve $e^{\hat\omega\theta} = R$ (Prop 2.5). Then solve $Av = p$ with $A = (I - e^{\hat\omega\theta})\hat\omega + \omega\omega^T\theta$. A is invertible for $\theta\neq0$, because its two pieces have mutually orthogonal null spaces.

**exp: se(3) → SE(3) is surjective, not injective** (inherited from SO(3)).
**Exponential coordinates of g:** $\xi\theta\in\mathbb R^6$.

### 5.6 Planar version, se(2) (from the review sheet) ★

$$
\hat\xi = \begin{bmatrix}\hat\omega & v\cr 0&0\end{bmatrix}\in\mathbb R^{3\times3},\quad \hat\omega = \begin{bmatrix}0&-\omega\cr \omega&0\end{bmatrix}\in so(2),\quad v\in\mathbb R^2,\quad \xi = (v,\omega)\in\mathbb R^3
$$

A rotation about the planar point $q = (q_x,q_y)$ has $v = -\hat\omega q = (\omega q_y, -\omega q_x)$.

**Every planar rigid motion is a pure rotation** about some point (possibly infinitely far away, which is the pure‑translation case).

---

## 6. Screws (L8) ★★

**Q:** for an arbitrary twist $(v,\omega)$, can we always find q with $v = -\omega\times q$?
**A: No.** $-\omega\times q$ is always ⟂ ω. If v has a component along ω, you need a **pitch** term:

$$
v = -\omega\times q + h\omega,\qquad \boxed{h = \frac{\omega^Tv}{\Vert\omega\Vert^2}}\ \ (\text{translation per radian})
$$

```
   Screw motion = rotate θ about the axis  +  slide hθ along the same axis

        q ●━━━━━━━━━━━━━━━━━━━━━━━━━━▶ ω  (screw axis)
             ↻ θ          ──hθ──▶

   h = 0 : pure rotation      h = ∞ : pure translation (ω = 0)
```

**Screw ↔ twist conversion:**

| Screw attribute | From the twist $\xi = (v,\omega)$ |
|---|---|
| Pitch | $h = \omega^Tv/\Vert\omega\Vert^2$ ($\infty$ if $\omega=0$) |
| Axis | $\lbrace\frac{\omega\times v}{\Vert\omega\Vert^2} + \lambda\omega\rbrace$ if $\omega\neq0$; $\lbrace0 + \lambda v\rbrace$ if $\omega=0$ |
| Magnitude | $M = \Vert\omega\Vert$ if $\omega\neq0$, else $\Vert v\Vert$ |

The point on the axis closest to the origin is $\omega\times v/\Vert\omega\Vert^2$.

**Unit twist and magnitude.** If $\xi_{unit}$ has $\Vert\omega\Vert= 1$ (or $\Vert v\Vert= 1$), then for $\xi = \xi_{unit}\theta$ the magnitude is $M = \Vert\omega\theta\Vert= |\theta|$. So **θ plays the role of the screw magnitude** (its absolute value, if you want to be picky about signs).

**Chasles (Thm 2.11):** *every* rigid motion in 3D is a screw motion. Instantaneous velocities are twists too, so **every velocity $\hat V^s$ or $\hat V^b$ has a screw** with zero, finite, or infinite pitch.

> **Worked example (HW2 P1 style, a different bolt).** A nut threads onto a bolt parallel to z through $(2, 0, 0)$ with 8 threads per inch.
> - One turn ($2\pi$ rad) moves it $1/8$ inch, so $h = \frac{1}{8\cdot2\pi} = \frac{1}{16\pi}$ in/rad.
> - Sign: for a right‑hand thread turned the right‑hand way about $+z$, it advances along $+z$.
>
> $$
> \xi = \begin{bmatrix}-z\times q + hz\cr  z\end{bmatrix},\quad -z\times(2,0,0) = -(0,2,0)\ \Rightarrow\ \xi = \Big(0,\ -2,\ \tfrac{1}{16\pi},\ 0,0,1\Big)
> $$

---

## 7. Forward kinematics with the Product of Exponentials (L8–L10) ★★★ (Q6)

### 7.1 The setup

- **Joint space** $Q$: one entry per joint ($S^1$ for a revolute joint, $\mathbb R$ for a prismatic one). Example: SCARA is $S^1\times S^1\times S^1\times\mathbb R = T^3\times\mathbb R$.
- **Joint i connects link i−1 to link i.** S is attached to the base (link 0), and T to the last link.
- **The forward kinematics map** is $g_{st}: Q\to SE(3)$.

### 7.2 The four‑step POE recipe ★★★

1. **Reference configuration** $\theta = 0$: write $g_{st}(0) = \begin{bmatrix}R_{st}(0) & p_{st}(0)\cr 0&1\end{bmatrix}$.
   - $p$ is the tool origin measured from S.
   - The columns of R are the tool axes in S coordinates.
   - Frames drawn parallel give $R = I$.
2. **Twist for each joint**, in S coordinates, at $\theta = 0$:
   - revolute: $\xi_i = (-\omega_i\times q_i, \omega_i)$
   - prismatic: $\xi_i = (v_i, 0)$
3. **Exponentials** $e^{\hat\xi_i\theta_i}$, using the shortcut from 5.3.
4. **Multiply:**

$$
\boxed{g_{st}(\theta) = e^{\hat\xi_1\theta_1}e^{\hat\xi_2\theta_2}\cdots e^{\hat\xi_n\theta_n}\,g_{st}(0)}
$$

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
> - $g_{st}(0) = \begin{bmatrix}I & (0,\ l_1+l_2,\ l_0)^T\cr 0&1\end{bmatrix}$
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

$$
e^{\hat\xi_1}e^{\hat\xi_2}e^{\hat\xi_3}\neq e^{\hat\xi_1+\hat\xi_2+\hat\xi_3}
$$

in general. $e^Ae^B = e^{A+B}$ only if A and B commute, and twists usually don't. (SE(3) is not abelian.)

**Rigid motions preserve distance.** $\Vert g\bar x\Vert= \Vert x\Vert$ for a vector x, and the difference of two points is a vector. So for example:

$$
\big\Vert e^{\hat\xi_1\theta_1}e^{\hat\xi_2\theta_2}\big(e^{\hat\xi_3\theta_3}\bar p - \bar q\big)\big\Vert= \big\Vert e^{\hat\xi_3\theta_3}\bar p - \bar q\big\Vert
$$

This is exactly the trick that drives SP3 in inverse kinematics.

**$e^{\hat\xi\theta}g = g$ is generally false.** Exponentials don't leave frames alone. They fix only points on their own axis (for a revolute joint).

**A tool point vs. a tool frame:**
- A point $p_t$ given in **tool coordinates** maps to its current spatial position by $g_{st}(\theta)\bar p_t$.
- A point given in **spatial coordinates at the reference configuration** moves by $e^{\hat\xi_1\theta_1}\cdots e^{\hat\xi_n\theta_n}\bar p_s(0)$, *not* by $g_{st}(\theta)\bar p_s(0)$. Watch which frame the point is written in.

---

## 8. The Adjoint transformation and order of exponentials (L11) ★★★

### 8.1 Definition

To change the coordinates a twist is written in, use

$$
\boxed{\mathrm{Ad}_g = \begin{bmatrix}R & \hat pR\cr 0&R\end{bmatrix}\in\mathbb R^{6\times6},\qquad \mathrm{Ad}_g^{-1} = \mathrm{Ad}_{g^{-1}} = \begin{bmatrix}R^T & -R^T\hat p\cr 0&R^T\end{bmatrix}}
$$

for $g = (p, R)$.

**Lemma 2.13:** $g\hat\xi g^{-1}$ is a twist, with coordinates $\mathrm{Ad}_g\xi$. Combined with 5.4:

$$
e^{(\mathrm{Ad}_g\xi)^\wedge\theta} = g\,e^{\hat\xi\theta}g^{-1}
$$

**What it does geometrically:** $\mathrm{Ad}_g\xi$ is the *same physical screw axis*, moved rigidly by g.

```
     original axis (ω, through q)          moved axis (Rω, through Rq + p)
            ↑ ω                                    ↗ ω' = Rω
            │                   ── g ──▶          ╱
       q ●──┘                               q' ●╱      q' = g(q) = Rq + p
```

### 8.2 The by‑hand shortcut (avoid 6×6 matrices) ★★★

For a **revolute** twist $\xi = (-\omega\times q, \omega)$:

$$
\boxed{\mathrm{Ad}_g\xi = \begin{bmatrix}-\omega'\times q'\cr \omega'\end{bmatrix},\quad \omega' = R\omega,\quad q' = Rq + p}
$$

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

$$
g_{st} = e^{\hat\xi_2'\theta_2}e^{\hat\xi_1\theta_1}g_{st}(0) = \big(e^{\hat\xi_1\theta_1}e^{\hat\xi_2\theta_2}e^{-\hat\xi_1\theta_1}\big)e^{\hat\xi_1\theta_1}g_{st}(0) = e^{\hat\xi_1\theta_1}e^{\hat\xi_2\theta_2}g_{st}(0)
$$

**Same answer.** The physical order doesn't matter, but if you move joint 1 first you have to use the *moved* twist.

**General rule:** the current spatial twist of joint i is

$$
\boxed{\xi_i' = \mathrm{Ad}_{e^{\hat\xi_1\theta_1}\cdots e^{\hat\xi_{i-1}\theta_{i-1}}}\,\xi_i}
$$

This depends only on the joints **before** i. Remember this; it is exactly column i of the spatial Jacobian (Section 13).

### 8.4 Expressing a twist in another frame (e.g., the tool frame) ★★★

If $\xi^s$ is written in S, then the same twist in frame B is $\xi^b = \mathrm{Ad}_{g_{sb}}^{-1}\xi^s$.

**The much easier way to do it:** a joint axis is **rigidly fixed to every link after it**. So the axis of joint i, written in the coordinates of any link frame $j\ge i$ (including the tool when j = n), is **constant**. You can read it straight off the reference configuration:

$$
\text{(joint } n\text{'s twist in tool coords, any } \theta) = \mathrm{Ad}^{-1}_{g_{st}(\theta)}\xi_n' = \mathrm{Ad}^{-1}_{g_{st}(0)}\xi_n
$$

**Algebraic proof:** $g_{st}(\theta) = (e^{1}\cdots e^{n-1})\,e^{n}\,g_{st}(0)$, and $e^{\hat\xi_n\theta_n}$ commutes with its own twist, so $\mathrm{Ad}_{e^{-\hat\xi_n\theta_n}}\xi_n = \xi_n$.

**In practice:** write ω and q of the last joint *as seen from the tool frame* in the picture, then form $(-\omega\times q, \omega)$. With frames parallel and the tool at distance $\ell$ along +y from the last axis point, for example, $q^T = (0, -\ell, 0)$.

---

## 9. Denavit–Hartenberg as relative twists (L9, L11) ★★★ (Q6b–d)

**This course's D‑H is not the textbook "a, α, d, θ" table.** It is POE written link to link:

$$
g_{st}(\theta) = g_{sl_1}(\theta_1)\,g_{l_1l_2}(\theta_2)\cdots g_{l_{n-1}l_n}(\theta_n)\,g_{l_nt},\qquad \boxed{g_{l_{i-1}l_i}(\theta_i) = e^{\hat\xi_{i-1,i}\theta_i}\,g_{l_{i-1}l_i}(0)}
$$

Here **$\xi_{i-1,i}$ is joint i's twist written in the coordinates of link frame $L_{i-1}$** (not S).
- Link frames sit at each link's proximal joint ($L_i$ on joint i's axis).
- $L_0 = S$, so $\xi_{0,1} = \xi_1$.

**Why it's constant:** axis i is fixed to link i−1, so its coordinates in $L_{i-1}$ never change, whatever the other joints do. You can define it at any configuration and get the same numbers.

**The relation to POE (L11):**

$$
\boxed{\xi_i = \mathrm{Ad}_{g_{sl_{i-1}}(0)}\,\xi_{i-1,i}\qquad\Longleftrightarrow\qquad \xi_{i-1,i} = \mathrm{Ad}^{-1}_{g_{sl_{i-1}}(0)}\,\xi_i}
$$

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

$$
g_{l_3l_4}(\theta_4) = \begin{bmatrix}R_x(\theta_4) & (I - R_x(\theta_4))q\cr 0&1\end{bmatrix}\begin{bmatrix}I & (0,l_1,0)^T\cr 0&1\end{bmatrix} = \begin{bmatrix}R_x(\theta_4) & (0,l_1,0)^T\cr 0&1\end{bmatrix},\quad q = (0,l_1,0)
$$

> **Shortcut ★★:** when $L_i$ sits **on** joint i's axis, the result is always
>
> $$
> g_{l_{i-1}l_i}(\theta_i) = \begin{bmatrix}e^{\hat\omega_i\theta_i} & p_{l_{i-1}l_i}(0)\cr 0&1\end{bmatrix}
> $$
>
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

$x = l_1c_1 + l_2c_{12}$, $y = l_1s_1 + l_2s_{12}$, $r = \sqrt{x^2+y^2}$, $\phi = \mathrm{atan2}(y,x)$.

$$
\cos\theta_2 = \frac{r^2 - l_1^2 - l_2^2}{2l_1l_2}\ \Rightarrow\ \theta_2 = \pm\cos^{-1}(\cdot)\qquad \theta_1 = \phi - \mathrm{atan2}(l_2s_2,\ l_1 + l_2c_2)
$$

That's 2 solutions (elbow up/down), 1 at full extension or full fold, and 0 if $r \gt l_1+l_2$ or $r \lt \lvert l_1-l_2\rvert$.
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
3. **Solvable iff** $\omega^Tu = \omega^Tv$ (same height) **and** $\Vert u'\Vert= \Vert v'\Vert$ (same radius)
4.

$$
\boxed{\theta = \mathrm{atan2}\big(\omega^T(u'\times v'),\ u'^Tv'\big)}
$$

5. **One solution**, or **infinitely many** if $u' = 0$ (p on the axis).

> **T/F ★★:** "A solution to $e^{\hat\xi\theta}p = q$ always exists as long as $p\neq q$." **False.** You need equal height and equal radius. Most pairs of points can't be reached by one rotation.

**By hand with a coordinate axis:** for $\omega = z$ and $r = 0$, θ is just $\mathrm{atan2}(q_y, q_x) - \mathrm{atan2}(p_y, p_x)$.

### 10.5 Subproblem 2: two intersecting axes ★★★

**Solve** $e^{\hat\xi_1\theta_1}e^{\hat\xi_2\theta_2}\bar p = \bar q$. The axes meet at r.

**Idea:** find the intermediate point $c = e^{\hat\xi_2\theta_2}p = e^{-\hat\xi_1\theta_1}q$. It must lie on **both**:
- the circle p sweeps around axis 2
- the circle q sweeps around axis 1

Two circles on a sphere of radius $\Vert p - r\Vert$ intersect in 0, 1, or 2 points.

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

$$
\alpha = \frac{(\omega_1^T\omega_2)\,\omega_2^Tu - \omega_1^Tv}{(\omega_1^T\omega_2)^2 - 1},\qquad \beta = \frac{(\omega_1^T\omega_2)\,\omega_1^Tv - \omega_2^Tu}{(\omega_1^T\omega_2)^2-1}
$$

$$
\gamma^2 = \frac{\Vert u\Vert^2 - \alpha^2 - \beta^2 - 2\alpha\beta\,\omega_1^T\omega_2}{\Vert\omega_1\times\omega_2\Vert^2}
$$

- $\gamma^2 \gt 0$ gives two solutions, $\gamma^2 = 0$ one, and $\gamma^2 \lt 0$ none.
- **Then do SP1 twice:** rotate p → c about $\omega_2$ to get $\theta_2$, and rotate c → q about $\omega_1$ to get $\theta_1$.

**Perpendicular axes** ($\omega_1^T\omega_2 = 0$), the common exam case:

$$
\alpha = \omega_1^Tv,\qquad \beta = \omega_2^Tu,\qquad \gamma^2 = \Vert u\Vert^2 - \alpha^2 - \beta^2
$$

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

**Solve** $\Vert\bar q - e^{\hat\xi\theta}\bar p\Vert= \delta$.

The circle of p around the axis meets the sphere of radius δ around q. That gives **0, 1, or 2** solutions.

**Steps:**
1. $u, v, u', v'$ as in SP1.
2. $\delta'^2 = \delta^2 - |\omega^T(p - q)|^2$ (remove the height difference)
3. $\theta_0 = \mathrm{atan2}(\omega^T(u'\times v'),\ u'^Tv')$
4.

$$
\boxed{\theta = \theta_0\pm\cos^{-1}\!\Big(\frac{\Vert u'\Vert^2 + \Vert v'\Vert^2 - \delta'^2}{2\Vert u'\Vert\Vert v'\Vert}\Big)}
$$

The last step is the law of cosines in the projected plane.

### 10.7 Full IK strategy: elbow manipulator (L12–13) ★★★

$e^1e^2e^3e^4e^5e^6 = g_1 := g_dg_{st}(0)^{-1}$, with wrist center $p_w$ (axes 4, 5, 6 meet there) and $p_b$ (axes 1, 2 meet there).

| Step | Trick | Subproblem | Number of solutions |
|---|---|---|---|
| **θ₃** | Apply to $p_w$: the wrist exponentials vanish, giving $e^1e^2e^3p_w = g_1p_w$. Subtract $p_b$ and take norms: $\Vert e^3p_w - p_b\Vert= \Vert g_1p_w - p_b\Vert= \delta$ | SP3 | 2 |
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

$$
\hat\omega^s_{ab} = \dot R_{ab}R_{ab}^T,\qquad \hat\omega^b_{ab} = R_{ab}^T\dot R_{ab},\qquad \omega^s = R_{ab}\,\omega^b
$$

> **T/F ★★★:** "If $R(t)\in SO(3)$, then $\dot R(t)\in so(3)$." **False.** $\dot R = \hat\omega^sR$, which is not skew in general. It's $\dot RR^T$ and $R^T\dot R$ that are skew. (The only exception is at $R = I$.)

**Velocity of a point on the body:** $v_{q_a} = \omega^s\times q_a$, or in body coordinates $v_{q_b} = \omega^b\times q_b$.

### 11.2 General (SE(3)) velocity ★★★

$$
\hat V^s_{ab} = \dot g_{ab}g_{ab}^{-1},\qquad V^s_{ab} = \begin{bmatrix}v^s\cr \omega^s\end{bmatrix} = \begin{bmatrix}-\dot RR^Tp + \dot p\cr  (\dot RR^T)^\vee\end{bmatrix}
$$

$$
\hat V^b_{ab} = g_{ab}^{-1}\dot g_{ab},\qquad V^b_{ab} = \begin{bmatrix}v^b\cr \omega^b\end{bmatrix} = \begin{bmatrix}R^T\dot p\cr  (R^T\dot R)^\vee\end{bmatrix}
$$

$$
\boxed{V^s = \mathrm{Ad}_{g_{ab}}V^b}\qquad V^b_{ab} = -V^s_{ba}\ \ (\text{Lemma 2.16})
$$

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

$$
\bar v_{q_a} = \hat V^s\bar q_a = \begin{bmatrix}\omega^s\times q_a + v^s\cr 0\end{bmatrix},\qquad \bar v_{q_b} = \hat V^b\bar q_b = \begin{bmatrix}\omega^b\times q_b + v^b\cr 0\end{bmatrix}
$$

### 11.3 The fastest way to get $V^s$ for a known motion ★★★

If the body **rotates about a fixed axis** (ω, through q) at rate $\dot\theta$, then $g(t) = e^{\hat\xi\theta(t)}g(0)$, and

$$
\boxed{V^s = \xi\,\dot\theta = \begin{bmatrix}-\omega\times q\cr \omega\end{bmatrix}\dot\theta}
$$

**The spatial velocity is just the joint twist times the joint rate.** Add translation along the axis to get a screw: $v^s = -\omega\times q\,\dot\theta + v_{\text{trans}}$.

**Integrating a constant velocity:**
- constant $V^s$: $g(t) = e^{\hat V^st}g(0)$
- constant $V^b$: $g(t) = g(0)e^{\hat V^bt}$

> **Worked example (HW3‑style, different numbers).** Body B starts with its origin at $(2,0,0)$ in A, axes parallel. It rotates about its own z‑axis at $\pi$ rad/s.
> - **Spatial velocity:** $\omega = (0,0,\pi)$, axis through $q = (2,0,0)$. $v^s = -\omega\times q = -\pi\,(z\times(2,0,0)) = -\pi(0,2,0)$, so $V^s = (0,-2\pi,0,\ 0,0,\pi)$.
> - **Point** $q_b = (1,0,0)$, so $q_a = (3,0,0)$. $v_{q_a} = \omega\times q_a + v^s = (0,3\pi,0) + (0,-2\pi,0) = (0,\pi,0)$.
>   - **Sanity check:** it's 1 unit from the spin axis and spinning at π, so its speed is π, tangent in +y ✓.
> - **Position after 0.5 s:** the body has rotated π/2 about the line x = 2, so $q_a = (2,0,0) + R_z(\tfrac\pi2)(1,0,0) = (2,1,0)$.
> - **Add translation at 3 units/s along +z:** $V^s = (0,-2\pi,3,\ 0,0,\pi)$, with pitch $h = \omega^Tv/\Vert\omega\Vert^2 = 3\pi/\pi^2 = 3/\pi$. After 0.5 s, $q_a = (2,1,1.5)$ (verified with $e^{\hat V^s\cdot0.5}$).

---

## 12. Wrenches (L14) ★★

**A wrench is a force and a moment at a point**, written in some frame's coordinates:

$$
F = \begin{bmatrix}f\cr \tau\end{bmatrix}\in\mathbb R^6
$$

For example, $F_b$ is applied at B's origin, in B's coordinates.

**Power pairs twists with wrenches:**

$$
\delta W = V^b\cdot F_b = v^b\cdot f_b + \omega^b\cdot\tau_b
$$

Two wrenches are **equivalent** if they do the same work for every motion.

**Changing frames uses the transpose of the Adjoint**, not a 4×4 g:

$$
\boxed{F_c = \mathrm{Ad}^T_{g_{bc}}F_b,\qquad \mathrm{Ad}_g^T = \begin{bmatrix}R^T & 0\cr  -R^T\hat p & R^T\end{bmatrix}\ \Rightarrow\ f_c = R^Tf_b,\quad \tau_c = R^T(\tau_b - p\times f_b)}
$$

Here p is the position of C's origin in B's coordinates.

**Where this comes from:** power is invariant. $V_c = \mathrm{Ad}^{-1}_{g_{bc}}V_b$, so $V_c^TF_c = V_b^TF_b$ forces $F_c = \mathrm{Ad}^T_{g_{bc}}F_b$.

**Physical meaning:** the force is the same vector (just rotated into new coordinates). The **moment changes** with the reference point, by $-p\times f$. So the **net wrench depends on where you measure it**: moving the reference point P changes the moment part.

> **T/F ★★:** "To move a wrench from frame A to frame B, multiply by a homogeneous transform in SE(3)." **False.** You need a 6×6 $\mathrm{Ad}^T$. (A 4×4 g acts on points, not wrenches.)

---

## 13. The manipulator Jacobian (MLS §3.4, *probably upcoming in lecture*; tested on both exams) ★★★ (Q7)

### 13.1 The spatial Jacobian

$$
V^s_{st} = J^s_{st}(\theta)\,\dot\theta,\qquad \boxed{J^s_{st}(\theta) = \big[\ \xi_1\ \ \xi_2'\ \ \cdots\ \ \xi_n'\ \big],\quad \xi_i' = \mathrm{Ad}_{e^{\hat\xi_1\theta_1}\cdots e^{\hat\xi_{i-1}\theta_{i-1}}}\,\xi_i}
$$

**Column i is joint i's twist, at the current configuration, in S coordinates.** You compute it with the shortcut from 8.2: move the axis point and direction using the earlier joints, then rebuild $(-\omega'\times q', \omega')$.
- Column 1 is always the reference twist $\xi_1$.
- **$J^s$ never depends on $\theta_n$.** Column i depends only on $\theta_1,\dots,\theta_{i-1}$.
- **At $\theta = 0$:** $J^s(0) = [\xi_1\ \cdots\ \xi_n]$, just the reference twists. That's why "assume θ = 0" makes exam problems easy.

### 13.2 The body Jacobian

$$
V^b_{st} = J^b_{st}\dot\theta,\qquad J^b_{st} = \big[\xi_1^\dagger\ \cdots\ \xi_n^\dagger\big],\quad \xi_i^\dagger = \mathrm{Ad}^{-1}_{e^{\hat\xi_i\theta_i}\cdots e^{\hat\xi_n\theta_n}g_{st}(0)}\,\xi_i
$$

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

$$
J^s_{st}(\theta) = \begin{bmatrix}0&0&-l_1s_1c_2\cr 0&0&l_1c_1c_2\cr 0&0&l_1s_2\cr 0&c_1&c_1\cr 0&s_1&s_1\cr 1&0&0\end{bmatrix}\qquad(\text{checked by differentiating } g_{st}(\theta)\text{ symbolically})
$$

**Velocity of a tool point at θ = 0** ($\bar q^s = (q_x,q_y,q_z,1)$):
- $V^s = J^s(0)\dot\theta = (0,\ l_1\dot\theta_3,\ 0,\ \dot\theta_2+\dot\theta_3,\ 0,\ \dot\theta_1)$
- $v_q = \omega\times q + v = \big(-\dot\theta_1q_y,\ \ \dot\theta_1q_x - (\dot\theta_2+\dot\theta_3)q_z + l_1\dot\theta_3,\ \ (\dot\theta_2+\dot\theta_3)q_y\big)$

**Body Jacobian at θ = 0:** $\xi_i^\dagger = \mathrm{Ad}^{-1}_{g_{st}(0)}\xi_i$. That's just each axis re‑measured from the tool origin $(0,0,l_1+l_2)$:
- $\xi_1^\dagger = (0,0,0,\ 0,0,1)$
- $\xi_2^\dagger$: axis through $(0,0,-(l_1+l_2))$ in tool coordinates, so $\xi_2^\dagger = (0,-(l_1+l_2),0,\ 1,0,0)$
- $\xi_3^\dagger = (0,-l_2,0,\ 1,0,0)$

**Joint torques for a body wrench $F_b = (f_x,f_y,f_z,\tau_x,\tau_y,\tau_z)$ at the tool, θ = 0:**

$$
\tau = (J^b_{st})^TF_b = \big(\tau_z,\ \ \tau_x - (l_1+l_2)f_y,\ \ \tau_x - l_2f_y\big)
$$

**Sanity check:** a +y force at height L about the x‑axis gives moment $(0,0,L)\times(0,f_y,0) = (-Lf_y, 0, 0)$ ✓.

### 13.4 Statics: $\tau = J^TF$ ★★★

$$
\boxed{\tau = (J^b_{st})^TF_t\ \ (\text{tool wrench in tool coords})\qquad\text{or}\qquad \tau = (J^s_{st})^TF_s\ \ (\text{spatial wrench})}
$$

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

## 14. True/False deep dive: the concepts behind every exam statement ★★★

The T/F blocks are 30 of the 70 points (20 statements × 1.5). Each statement usually hinges on **one hidden assumption**: a unit‑length requirement, a "for all" vs. "for some", which frame something lives in, or which operation a group uses.

**How to attack any T/F statement:**

1. **Translate the words into math.** "A relative to C" becomes $R_{CA}$. "Rotate joint i" becomes "θᵢ changes, everything else fixed."
2. **Find the quantifier.** "Always / for all" means one counterexample makes it **False**. "There exists / for some" means one example makes it **True**.
3. **Check the fine print of the formula it quotes.** Does it need $\Vert\omega\Vert= 1$? $\omega \neq 0$? Pitch 0? Which frame are the coordinates in?
4. **Try the simplest case:** $R = I$, $\theta = 0$, an axis through the origin, a 1‑ or 2‑joint robot. Most counterexamples are trivial ones.
5. **Do a dimension check.** A 6×n matrix times a 6‑vector is undefined, so a formula like that is immediately False.

Each entry below gives the statement, the answer, **why**, and **how the exam might twist it**.

---

### 14.1 Rotations and SO(3)

**1. $R_{CA} = R_{BC}^T R_{BA}$.** — **True**

Chain the subscripts so inner letters match: $R_{CA} = R_{CB}R_{BA}$. Reversing the subscripts inverts the rotation, and the inverse of a rotation is its transpose, so $R_{CB} = R_{BC}^T$.

*Exam twist:* the statement is written in words, for example "frame A has orientation $R_1$ relative to B, frame C has orientation $R_2$ relative to B; then A relative to C is $R_2^TR_1$." The convention is **"X relative to Y" = $R_{YX}$** (the reference frame comes first). So $R_1 = R_{BA}$, $R_2 = R_{BC}$, and $R_{CA} = R_2^TR_1$. That's True. If they swapped to $R_1^TR_2$, you'd get $R_{AB}R_{BC} = R_{AC}$, which is C relative to A, so False.

**Method that never fails:** write the target subscripts, write the givens as subscripts, connect them, and put a transpose on every given you had to traverse backwards.

**2. so(3) is an abelian group under matrix addition.** — **True**

Check the axioms:
- **Closure:** $(A+B)^T = -A - B$, so the sum of two skew matrices is skew.
- **Identity:** the zero matrix is skew.
- **Inverse:** $-A$ is skew.
- **Associativity:** inherited from matrix addition.
- **Commutativity:** $A + B = B + A$ always.

In fact *any* vector space is an abelian group under addition, and so(3) is a 3‑dimensional vector space.

*Exam twist:* the statement names the operation. Under **addition**, True. Under **multiplication**, False (next item).

**3. so(3) is a group under matrix multiplication.** — **False**

It fails three axioms at once:
- **No closure.** Use the identity $\hat a\hat b = ba^T - (a^Tb)I$. Then $\hat x\hat y = yx^T$, a matrix with a single 1 in position (2,1). That isn't skew.
- **No identity.** $I$ is not skew ($I^T = I \neq -I$).
- **No inverses.** Every 3×3 skew matrix has determinant 0 (item 5), so none of them is invertible.

**4. SO(2) is abelian; SO(3) is not.** — **True**

- **SO(2):** every element is $R(\theta)$, and $R(a)R(b) = R(a+b) = R(b)R(a)$. Planar rotations all share the same axis (out of the plane), so they commute.
- **SO(3) counterexample:** apply both orders to the vector $x$.
  - $R_x(90°)R_z(90°)\,x$: first $R_z$ sends $x \to y$, then $R_x$ sends $y \to z$. Result: $z$.
  - $R_z(90°)R_x(90°)\,x$: first $R_x$ leaves $x$ alone, then $R_z$ sends $x \to y$. Result: $y$.

  Different results, so they don't commute.

*Physical version:* take a book, rotate it 90° about x then 90° about z, and compare with the other order.

*Related:* SE(2) and SE(3) are also non‑abelian. Even in the plane, rotating then translating ≠ translating then rotating.

**5. If $\Vert\omega\Vert= 1$, then $\hat\omega$ is invertible and $\det\hat\omega = 1$.** — **False**

$\hat\omega\,\omega = \omega\times\omega = 0$, so ω is a nonzero vector in the null space of $\hat\omega$. That makes $\hat\omega$ singular: $\det\hat\omega = 0$.

**Second proof, for any 3×3 skew matrix S:** $\det S = \det S^T = \det(-S) = (-1)^3\det S = -\det S$, so $\det S = 0$.

- $\hat\omega$ has **rank 2** when $\omega \neq 0$. Its eigenvalues are $0, \pm i\Vert\omega\Vert$, and their product is 0.

*Exam twist:* the unit‑length condition is a distractor. It matters for Rodrigues' formula but not here.

*Related:* this is why you can't "solve $\hat\omega x = b$" uniquely. You can only recover ω from $\hat\omega$ with the vee map.

**6. Rodrigues' formula $e^{\hat\omega\theta} = I + \hat\omega\sin\theta + \hat\omega^2(1-\cos\theta)$ holds for every $\hat\omega\in so(3)$ and every θ.** — **False**

The derivation used $\hat\omega^3 = -\hat\omega$ (Lemma 2.3 with $\Vert\omega\Vert= 1$) to collapse the series into sin and cos. For non‑unit ω, $\hat\omega^3 = -\Vert\omega\Vert^2\hat\omega$, so the coefficients come out wrong.

**Counterexample:** $\omega = (0,0,2)$, $\theta = \pi/2$.
- The true answer is a rotation by $2\cdot\frac{\pi}{2} = \pi$ about z: $\mathrm{diag}(-1,-1,1)$.
- The formula gives $I + \hat\omega + \hat\omega^2$, with $\hat\omega^2 = \mathrm{diag}(-4,-4,0)$:

$$
\begin{bmatrix}-3&-2&0\cr 2&-3&0\cr 0&0&1\end{bmatrix}
$$

That isn't even orthogonal.

**The fix:** normalize first. $\hat\omega\theta = \widehat{\big(\tfrac{\omega}{\Vert\omega\Vert}\big)}\,(\Vert\omega\Vert\theta)$.

*Exam twist:* the 2026 exam phrased it the other way: "Rodrigues' formula **does not hold for some** $\hat\omega\in so(3)$ and θ." That is **True**, by the same counterexample. Watch for negations.

**7. exp: so(3) → SO(3) is surjective.** — **True**

This is Prop 2.5, proved by construction. Given any R:
- $\theta = \cos^{-1}\frac{\mathrm{tr}R - 1}{2}$
- ω from the skew part (or the special formulas at θ = 0 or π)

*Physical meaning:* **Euler's theorem.** Every orientation can be reached by one rotation about one fixed axis. No matter how you tumble a body, the net result is a single axis‑angle rotation.

**8. exp: so(3) → SO(3) is injective.** — **False**

Many exponential coordinates give the same R:
- $(\omega,\theta)$
- $(-\omega,-\theta)$
- $(-\omega, 2\pi - \theta)$
- $(\omega, \theta + 2\pi k)$
- at $R = I$: $(\text{any }\omega, 0)$

**Picture:** rotating +90° about ω is the same as rotating +270° about −ω.

*Exam twist:* "exp is one‑to‑one", "exp is a bijection", or "exp has an inverse": all False. "exp is onto": True. The same pattern holds for se(3) → SE(3) and so(2) → SO(2).

**9. There exists R ∈ SO(3) with infinitely many exponential coordinates (ω, θ) with θ ∈ [0, 2π).** — **True**

Take $R = I$: θ = 0, and *every* unit ω satisfies $e^{\hat\omega\cdot 0} = I$.

Why the interval matters: restricting θ to $[0, 2\pi)$ kills the $2\pi k$ copies. After that:
- **$R \neq I$, θ ∈ (0, π):** exactly two solutions, $(\omega,\theta)$ and $(-\omega, 2\pi-\theta)$.
- **θ = π:** two solutions, $(\pm\omega, \pi)$.
- **Only $R = I$ gives infinitely many.**

So the statement is an "exists" question, and $R = I$ is the witness.

**10. θ = π is a singularity of the axis‑angle (exponential coordinate) representation.** — **False**

A **singularity** means the parameterization loses information: the parameters can't be recovered smoothly or uniquely from R. At θ = 0 that happens, because ω is completely undefined.

At θ = π, ω *is* determined (up to the harmless sign, since ±ω by π are the same rotation). Formula (2.18) only breaks because it divides by $\sin\pi = 0$. That's a breakdown of one *formula*, not of the representation (HW1 P3). Use instead:

$$
R = 2\omega\omega^T - I\ \Rightarrow\ \omega\omega^T = \tfrac{1}{2}(R+I)
$$

So $\omega_i = \pm\sqrt{(r_{ii}+1)/2}$, with the relative signs taken from $r_{ij} = 2\omega_i\omega_j$.

**11. If $A^{-1} = A^T$ and the eigenvalues of A are $1, e^{-i\pi/2}, e^{i\pi/2}$, then A ∈ SO(3).** — **True**

- Orthogonality gives $\det A = \pm 1$.
- The determinant is the product of the eigenvalues: $1\cdot e^{-i\pi/2}e^{i\pi/2} = 1$.
- Orthogonal with determinant +1 is exactly the definition of SO(3).

**Geometric reading:** eigenvalue 1 is the axis, and $e^{\pm i\pi/2}$ means a rotation by 90°.

*Exam twist:*
- Eigenvalues $\lbrace-1, e^{\pm i\phi}\rbrace$ give $\det = -1$: a **rotoreflection**, not in SO(3).
- If the statement drops orthogonality, eigenvalues alone aren't enough. A non‑orthogonal matrix can have those eigenvalues (any $P\,\mathrm{diag}(\cdot)P^{-1}$).

**12. If R(t) ∈ SO(3), then $\dot R(t)\in so(3)$.** — **False**

Differentiate $RR^T = I$: $\dot RR^T + R\dot R^T = 0$. So **$\dot RR^T$ is skew** (and likewise $R^T\dot R$), but $\dot R$ itself is $\dot R = \hat\omega^s R$: a skew matrix *times* a rotation. That's not skew in general.

**Counterexample:** $R(t) = R_z(t)$ at $t = \pi/2$:

$$
\dot R = \begin{bmatrix}-\sin t&-\cos t&0\cr \cos t&-\sin t&0\cr 0&0&0\end{bmatrix}_{t=\pi/2} = \begin{bmatrix}-1&0&0\cr 0&-1&0\cr 0&0&0\end{bmatrix}
$$

That's symmetric, not skew.

- $\dot R$ *is* skew at the instant $R = I$, which is why this is easy to believe.
- **Geometric picture:** so(3) is the tangent space at the **identity**. The tangent space at R is $so(3)\cdot R$.

**13. The map from Euler angles (α, β, γ) to R ∈ SO(3) is injective.** — **False**

Three separate reasons:
1. **Singularity:** at β = 0, $R = R_z(\alpha+\gamma)$, so only the sum matters. For $R = I$, every $(\alpha, 0, -\alpha)$ works.
2. **Periodicity:** adding 2π to any angle gives the same R.
3. **The twin solution:** $(\alpha+\pi, -\beta, \gamma+\pi)$ gives the same R as $(\alpha,\beta,\gamma)$. This is why the inverse uses $\mathrm{atan2}(+\sqrt{\cdot}, \cdot)$ and picks one branch.

**Topological fact (L4):** *no* 3‑parameter representation of SO(3) is both global and singularity‑free. So any "Euler angles are global / one‑to‑one / singularity‑free" statement is False.

**14. ZYX (yaw‑pitch‑roll) Euler angles are singular at pitch θ = ±π/2.** — **True**

At θ = π/2 the matrix becomes

$$
R = \begin{bmatrix}0 & \sin(\psi-\phi) & \cos(\psi-\phi)\cr  0 & \cos(\psi-\phi) & -\sin(\psi-\phi)\cr  -1 & 0 & 0\end{bmatrix}
$$

which depends **only on ψ − φ**. Yaw and roll have become rotations about the same physical axis, so you lose a degree of freedom. This is **gimbal lock** (an airplane pointing straight down).

*Compare:* ZYZ is singular at β = 0 or π instead. Know which convention has which singularity.

---

### 14.2 SE(3), twists, and screws

**15. The bottom‑right element of a twist $\hat\xi\in se(3)$ is 1.** — **False**

$$
\hat\xi = \begin{bmatrix}\hat\omega & v\cr 0&0\end{bmatrix}
$$

The **entire bottom row is zero**. Two ways to see it:
- **Twists are velocities.** $\hat\xi = \dot gg^{-1}$, and the bottom row of g is the constant $(0,0,0,1)$, whose derivative is 0.
- **se(3) is a vector space.** If the corner were 1, adding two twists would give 2, so the sum wouldn't be a twist.

*Don't confuse it with* $\bar g\in SE(3)$, whose corner **is** 1.

**16. For a twist $\hat\xi$ and homogeneous point $\bar p$, $\hat\xi\bar p$ is a homogeneous vector (bottom entry 0).** — **True**

$$
\hat\xi\bar p = \begin{bmatrix}\omega\times p + v\cr  0\end{bmatrix}
$$

That's the **velocity** of the point p under the motion ξ, and velocities are vectors. The bottom entry is $0\cdot p + 0\cdot 1 = 0$.

**17. For any twist ξ = (v, ω) with ω ≠ 0, v can be written as v = −ω × q for a point q on the axis.** — **False**

$-\omega\times q$ is **always perpendicular to ω**. If v has a component along ω, no q can produce it. Decompose instead:

$$
v = \underbrace{v - (\omega^Tv)\omega}_{=\,-\omega\times q,\ \ q = \omega\times v} + \underbrace{(\omega^Tv)\,\omega}_{h\omega}\qquad(\Vert\omega\Vert= 1)
$$

The general twist is a **screw**: $v = -\omega\times q + h\omega$ with pitch $h = \omega^Tv$.

**Counterexample:** $\xi = (0,0,1,\ 0,0,1)$. Here v is parallel to ω, with pitch 1: rotating about z while moving up z.

*Exam twist:* "for a **revolute joint**, v = −ω × q": True. That's the pitch‑zero case.

**18. For a revolute joint whose axis passes through the spatial origin, v = 0.** — **True**

Choose q = 0 on the axis: $v = -\omega\times 0 = 0$. Any other q on the axis gives the same v (v doesn't depend on which axis point you pick).

- The **converse** also holds: a revolute twist with v = 0 has its axis through the origin.
- This **doesn't** apply to prismatic joints. Their twist is $(v,0)$ with $\Vert v\Vert= 1$, so v is never 0.

**19. A unit twist ξ = (v, ω) has 2‑norm 1, for pure rotation and pure translation alike.** — **False**

"Unit" means **$\Vert\omega\Vert= 1$** (or $\Vert v\Vert= 1$ when ω = 0). It says nothing about the whole 6‑vector.
- **Pure translation:** $(v, 0)$ with $\Vert v\Vert= 1$, so the norm is 1. ✓
- **Pure rotation:** a rotation about z through $q = (0,2,0)$ has $v = q\times z = (2,0,0)$, so $\Vert\xi\Vert= \sqrt{4+1} = \sqrt5$. ✗

The statement requires both, so it's False. The 2‑norm only equals 1 when the rotation axis passes through the origin.

*Deeper reason:* v carries units of length × (1/time) and ω carries 1/time, so the Euclidean norm of a twist isn't even frame‑invariant.

**20. The formula $e^{\hat\xi\theta} = \begin{bmatrix}e^{\hat\omega\theta} & (I-e^{\hat\omega\theta})(\omega\times v) + \omega\omega^Tv\theta\cr 0&1\end{bmatrix}$ applies to any twist, including pure translation.** — **False**

Plug in ω = 0: $e^{0} = I$ and the top‑right becomes $(I - I)(0\times v) + 0 = 0$. That gives **the identity**, but the correct answer is

$$
\begin{bmatrix}I & v\theta\cr 0&1\end{bmatrix}
$$

The formula's derivation also assumed $\Vert\omega\Vert= 1$ (shifting to the axis and normalizing). Pure translation has no axis to shift to, so it gets **its own case**.

**21. If g ∈ SE(3) is a pure translation, then det g = 0.** — **False**

g is block upper‑triangular, so $\det g = \det R\cdot 1 = \det I = 1$. Every element of SE(3) is invertible, with $g^{-1} = (−R^Tp,\ R^T)$, so its determinant can never be 0.

*Where the confusion comes from:* the **twist** $\hat\xi$ always has determinant 0 (zero bottom row), and $\hat\omega$ does too.

**22. The exponential map (ξ̂, θ) ↦ g ∈ SE(3) is surjective.** — **True**

This is Prop 2.9, proved by construction:
- If R = I: $v = p/\Vert p\Vert$, $\theta = \Vert p\Vert$.
- Otherwise: get (ω, θ) from R, then solve $Av = p$ (A is invertible for θ ≠ 0).

*Physical meaning:* **Chasles' theorem.** Every rigid displacement is a single screw motion.

*Exam twist:* "injective" is False, for the same reasons as SO(3).

**23. If M is the magnitude of the screw for twist $\hat\xi\theta$ with ξ a unit twist, then M = θ.** — **True** (in the course's convention)

- Magnitude is $M = \Vert\omega\Vert$ for the twist $\xi\theta$ (or $\Vert v\Vert$ if ω = 0).
- With ξ unit, the ω‑part of $\xi\theta$ is $\omega_u\theta$, so $M = \Vert\omega_u\theta\Vert= |\theta|$.
- Lecture 8 states that "θ is the magnitude of ξ".

**Pedantic caveat:** a magnitude can't be negative, so strictly $M = |\theta|$. Go with the lecture's definition (**True**) unless the statement deliberately hinges on θ < 0.

**24. $\exp(\hat\xi_1)\exp(\hat\xi_2)\exp(\hat\xi_3) = \exp(\hat\xi_1 + \hat\xi_2 + \hat\xi_3)$.** — **False**

$e^Ae^B = e^{A+B}$ holds when A and B **commute**. Twists generally don't (SE(3) isn't abelian).

**Concrete check with rotations:** $R_x(90°)R_z(90°)$ has trace 0, so it's a rotation by **120°**. But $e^{(\hat x+\hat z)\pi/2}$ is a rotation about $(x+z)/\sqrt2$ by $\frac{\pi}{2}\sqrt2 \approx$ **127°**. Different.

**When it *does* hold:** twists that commute, such as a revolute and a prismatic joint on the same line, or two rotations about the same axis.

---

### 14.3 Forward kinematics

**25. To use $g_{st}(\theta) = e^{\hat\xi_1\theta_1}\cdots e^{\hat\xi_n\theta_n}g_{st}(0)$, the robot must physically move its joints in reverse order.** — **False**

**The final pose depends only on the joint values, not on the order the motors moved.** Reading the product right to left as "move joint n first" is just a *bookkeeping* device. It lets every twist stay at its reference value, because joints 1…n−1 haven't moved yet when joint n moves.

If you instead imagine moving joint 1 first, you must use the moved axes, $\xi_2' = \mathrm{Ad}_{e^{\hat\xi_1\theta_1}}\xi_2$, and so on. L11 shows the product collapses to the same formula.

**26. For a 3‑DOF robot, $\xi_3' = \mathrm{Ad}_{e^{\hat\xi_1\theta_1}e^{\hat\xi_2\theta_2}}\xi_3$ is the third joint's twist at the current configuration, in the spatial frame.** — **True**

Joints 1 and 2 carry axis 3 rigidly, by $g = e^{\hat\xi_1\theta_1}e^{\hat\xi_2\theta_2}$, and the Adjoint moves a twist by g.
- Joint 3's own rotation doesn't move axis 3, so $e^{\hat\xi_3\theta_3}$ is correctly **absent**.
- Joints after 3 don't affect it either.

This is exactly **column 3 of the spatial Jacobian**.

*Exam twist:* a version with $e^{\hat\xi_3\theta_3}$ included is still numerically correct (a twist is invariant under its own exponential). A version with the exponentials in the wrong order, or one that includes later joints, is False.

**27. If g ∈ SE(3) is the configuration of joint i's link frame, then $\exp(\hat\xi_i\theta_i)g = g$.** — **False**

A revolute exponential fixes only the **points on its axis**. A frame is a point (its origin) *plus* three directions. Even if the origin sits on the axis, the frame's x, y, z directions rotate.
- So $e^{\hat\xi_i\theta_i}g = g$ only when θᵢ = 0 (mod 2π).
- A prismatic exponential fixes nothing unless θ = 0.

*What is true:* $e^{\hat\xi\theta}\bar p = \bar p$ for points p **on** the axis.

**28. Let p be a point attached to the tool, given in spatial coordinates at the reference configuration. Then $g_{st}(\theta)\bar p$ is its current location.** — **False**

$g_{st}(\theta)$ maps **tool coordinates** to spatial coordinates. For a point given in **spatial** coordinates at θ = 0, the right formula is

$$
\bar p(\theta) = e^{\hat\xi_1\theta_1}\cdots e^{\hat\xi_n\theta_n}\,\bar p(0) = g_{st}(\theta)\,g_{st}(0)^{-1}\,\bar p(0)
$$

**Quick counterexample:** at θ = 0, take p at the tool origin, $(0, L, 0)$, with $g_{st}(0) = [I \mid (0,L,0)]$. Then $g_{st}(0)\bar p = (0, 2L, 0)$, double‑counted. That's wrong even before anything moves.

*What would be True:* "p given in **tool** coordinates; $g_{st}(\theta)\bar p_t$ gives its spatial location."

**29. For homogeneous points p, q and twists ξ₁, ξ₂, ξ₃: $\Vert e^{\hat\xi_1\theta_1}e^{\hat\xi_2\theta_2}(e^{\hat\xi_3\theta_3}\bar p - \bar q)\Vert= \Vert e^{\hat\xi_3\theta_3}\bar p - \bar q\Vert$.** — **True**

- Point minus point is a **vector** (bottom entry 0).
- A rigid motion acts on a vector by rotation only: $g\bar x = (Rx, 0)$.
- Rotations preserve length.

**This is the engine of Subproblem 3**: strip off the leading exponentials by taking norms.

*Trap version:* $\Vert g\bar p\Vert= \Vert\bar p\Vert$ for a single **point** is False. Translation changes a point's distance from the origin.

**30. The relative transformation between link frames L₀ and L₁ is $g_{l_0l_1}(\theta) = e^{\hat\xi_1\theta_1}g_{l_0l_1}(0)$, where $\hat\xi_1$ is the standard POE twist for joint 1.** — **True**

The D‑H formula is $g_{l_{i-1}l_i} = e^{\hat\xi_{i-1,i}\theta_i}g_{l_{i-1}l_i}(0)$, with the twist in $L_{i-1}$ coordinates. For i = 1, $L_0 = S$, so $\xi_{0,1} = \xi_1$, the POE twist.

*Exam twist:* the same statement for **i = 2 or higher** using the POE twist $\xi_i$ is False in general. You need $\xi_{i-1,i} = \mathrm{Ad}^{-1}_{g_{sl_{i-1}}(0)}\xi_i$ (unless $L_{i-1}$ coincides with S).

**31. The D‑H twist $\xi_{i-1,i}$ (joint i's twist in $L_{i-1}$ coordinates) can be defined at the current, non‑zero configuration.** — **True**

Joint i's axis is rigidly attached to link i−1, and so is frame $L_{i-1}$. Their relative geometry never changes, so $\xi_{i-1,i}$ is the **same at every configuration**:
- joints 1…i−1 move both together
- θᵢ rotates about the axis itself

You can read it off whichever configuration is easiest.

*Contrast:* the **POE** twist $\xi_i$ must be read at θ = 0. Its current version $\xi_i'$ changes with θ₁…θᵢ₋₁.

---

### 14.4 Inverse kinematics

**32. Given a robot and $g_d$, the maximum number of IK solutions is independent of how the IK is solved analytically.** — **True**

The solution set $\lbrace\theta : g_{st}(\theta) = g_d\rbrace$ is fixed by the robot's kinematics and $g_d$. A method can **miss** branches (a bad derivation might find 4 of 8), but it can't change how many exist.

*Exam twist:* "different choices of subproblems give different numbers of solutions": False.

**33. For a planar RRRR manipulator and $g_d\in SE(2)$, the maximum number of IK solutions is 4.** — **False**

SE(2) has dimension 3 ($x, y, \phi$), but the robot has 4 joints. Four unknowns and three equations leave a **1‑parameter family** of solutions: **infinitely many**. This is a "self‑motion": the elbow can slide around while the tool stays still.

| Robot | Solutions |
|---|---|
| Planar RR (position only) | at most 2 |
| Planar RRR in SE(2) | at most 2 |
| Planar RRRR in SE(2) | ∞ (redundant) |
| 7‑DOF arm in SE(3) (e.g., Canadarm2) | ∞ (redundant) |

**34. If ξ is a zero‑pitch twist and p lies on its axis, then $e^{\hat\xi\theta}\bar p = \bar p$ for all θ.** — **True**

Use the shortcut form: $e^{\hat\xi\theta}\bar p = R(p - q) + q$. If p is on the axis, $p - q$ is parallel to ω, and R leaves ω‑direction vectors unchanged. So the result is $p - q + q = p$.

This is **the** IK trick for killing exponentials.

*Fails if:*
- the pitch is nonzero (the screw slides p along the axis)
- the joint is prismatic (everything translates)

**35. Given points p, q and a zero‑pitch twist ξ, a solution θ to $e^{\hat\xi\theta}\bar p = \bar q$ always exists as long as p ≠ q.** — **False**

A rotation keeps p on its circle around the axis. So q must be:
1. at the **same height** along ω: $\omega^T(p - r) = \omega^T(q - r)$
2. at the **same radius** from the axis: $\Vert u'\Vert= \Vert v'\Vert$

**Counterexamples** (ω = z through the origin, p = (1,0,0)):
- q = (0,2,0): radius 1 vs. 2, so no solution.
- q = (0,1,5): height 0 vs. 5, so no solution.

The "p ≠ q" condition is a red herring. When p = q, θ = 0 works trivially.

**36. A manipulator at the outer boundary of its workspace (fully extended) has only a single IK solution.** — **False** (in general)

- **Planar RR:** at full extension the elbow‑up and elbow‑down branches merge into one. That's where the intuition comes from.
- **Spatial arms:** other branches survive.
  - **Shoulder flip:** rotate the base by π and swing the shoulder "over the top" to reach the same point.
  - **Wrist flip:** $(\theta_4+\pi, -\theta_5, \theta_6+\pi)$ gives the same tool orientation.

So the count drops, but not to one.

**37. Swapping the order of a revolute joint and a prismatic joint about the same axis (e.g., joints 3 and 4) leaves the IK solutions unchanged.** — **True**

Rotation about a line and translation along the same line **commute**. With $\hat\xi_R = \begin{bmatrix}\hat\omega & -\omega\times q\cr 0&0\end{bmatrix}$ and $\hat\xi_P = \begin{bmatrix}0&\omega\cr 0&0\end{bmatrix}$, both products $\hat\xi_R\hat\xi_P$ and $\hat\xi_P\hat\xi_R$ are 0, so $e^{\hat\xi_R\theta}e^{\hat\xi_Pd} = e^{\hat\xi_Pd}e^{\hat\xi_R\theta}$.

So $g_{st}$ is the same function of (θ, d), and the solutions are the same.

*Physical picture:* a screwdriver. Push then twist, or twist then push, and you end up in the same place.

**Related: is a single revolute angle unique given all the other joints?** — **Yes**, in $[0, 2\pi)$. A rotation angle about a fixed axis is determined uniquely by where it sends any point *off* the axis (SP1 has one solution). So "θ₃ has multiple solutions given the other joints" is False.

---

### 14.5 Velocities, wrenches, Jacobians

**38. Every instantaneous velocity $\hat V^s_{st}$ or $\hat V^b_{st}$ has an associated screw with pitch zero, finite, or infinite.** — **True**

$V^s$ and $V^b$ are twists. Every twist (v, ω) defines:
- an axis
- a pitch $h = \omega^Tv/\Vert\omega\Vert^2$ (∞ if ω = 0)
- a magnitude

Physically, at each instant the body is screwing about an **instantaneous screw axis**:

| Pitch | Motion |
|---|---|
| 0 | pure spin |
| ∞ | pure translation |
| finite, nonzero | helical |

**39. The linear part $v^s$ of the spatial velocity is the velocity of the body frame's origin.** — **False**

$v^s$ is the velocity of the (imaginary) **body point currently at the spatial origin A**. The velocity of B's origin p is

$$
\dot p = \omega^s\times p + v^s
$$

**Counterexample (Section 11's worked example):** a body spins about a z‑axis through (2,0,0), and its origin sits **on** that axis. So its origin isn't moving ($\dot p = 0$), yet $v^s = (0,-2\pi,0) \neq 0$.

*What is the origin's velocity?* $v^b$ (next item), or $\dot p$ in A coordinates.

**40. Since $v_{q_b} = \hat V^b_{ab}\bar q_b$, the velocity of the body frame's origin is not necessarily zero in body coordinates.** — **True**

Plug in $\bar q_b = (0,0,0,1)$: $\hat V^b\bar 0 = (v^b, 0)$, and $v^b = R^T\dot p$. That's the origin's velocity **relative to A**, just *expressed* in B's axes. It's zero only if the origin is stationary.

*The trap:* "B can't move relative to itself, so it's zero." But $V^b_{ab}$ is B's motion **relative to A**. Only the *coordinates* are B's.

**41. To express a wrench $F_a$ in frame B, multiply it by a homogeneous transformation in SE(3).** — **False**

A wrench is a 6‑vector $(f, \tau)$, and a 4×4 g acts on 4‑vectors (points and vectors), so the dimensions don't even match. More importantly, the **moment changes** when you move the reference point, $\tau_{new} = \tau - p\times f$. That mixes force into moment, which a 4×4 can't capture. The correct transform is 6×6:

$$
F_c = \mathrm{Ad}^T_{g_{bc}}F_b
$$

**42. If n wrenches are applied to a rigid body at frames $C_i$, the net wrench $F_p$ with respect to frame P doesn't depend on where P is on the body.** — **False**

- The net **force** is the same everywhere: $f_p = \sum R^Tf_i$.
- The net **moment** depends on the reference point: $\tau_P = \tau_O - p\times f$.

**Example:** a force (0,0,−10) applied at (1,0,0).
- About the origin, the moment is $(1,0,0)\times(0,0,-10) = (0,10,0)$.
- About (1,0,0), the moment is 0.

*Only exception:* a pure couple (net f = 0) has the same moment everywhere.

**43. The spatial Jacobian $J^s_{st}(\theta)$ does not depend on θₙ.** — **True**

Column i is $\mathrm{Ad}_{e^{\hat\xi_1\theta_1}\cdots e^{\hat\xi_{i-1}\theta_{i-1}}}\xi_i$, which uses only θ₁…θᵢ₋₁. Even the last column uses only up to θₙ₋₁.

**Physically:** the last joint spins the tool about axis n, but it doesn't move any joint axis in space.

**44. The body Jacobian $J^b_{st}(\theta)$ does not depend on θₙ.** — **False**

Column i is $\xi_i^\dagger = \mathrm{Ad}^{-1}_{e^{\hat\xi_{i+1}\theta_{i+1}}\cdots e^{\hat\xi_n\theta_n}g_{st}(0)}\xi_i$, which **does** contain θₙ for every i < n. The body Jacobian is instead independent of **θ₁**: rotating joint 1 carries the whole arm, tool included, rigidly, so nothing changes relative to the tool.

**Memory hook:**
- spatial ↔ independent of the **last** joint
- body ↔ independent of the **first** joint

*Exam twist:* "the manipulator Jacobian (spatial or body) doesn't depend on θₙ": False, because the body one does depend on it.

**45. Rotating the i‑th joint changes the i‑th column of the body Jacobian.** — **False**

$\xi_i^\dagger$ depends only on θᵢ₊₁…θₙ, since $\mathrm{Ad}_{e^{-\hat\xi_i\theta_i}}\xi_i = \xi_i$ (a twist is invariant under its own exponential).

**Physically:** rotating joint i spins the tool around axis i. Seen from the tool, axis i hasn't moved, because the tool orbited *around* it.

**Bonus:** the same holds for the spatial Jacobian. Column i of $J^s$ depends only on θ₁…θᵢ₋₁. **Moving joint i never changes column i of either Jacobian.**

**46. For a body wrench F on link 3 of a 5‑link manipulator, the joint torque to resist it is $\tau = J^b_{st}(\theta)F$.** — **False**

Two errors (either one is enough):
1. **The transpose is missing.** $J$ is 6×5, so $J F$ with a 6‑vector F is undefined. The statics law is $\tau = J^TF$.
2. **It's the wrong Jacobian.** A wrench on **link 3** needs link 3's body Jacobian $J^b_{sl_3}$. Its columns for joints 4 and 5 are **zero**, because joints after link 3 can't push on link 3.

The tool Jacobian $J^b_{st}$ belongs to wrenches at the tool. (There's also a sign: *resisting* an external wrench uses $-J^TF$.)

**47. At θ = 0, the spatial Jacobian is $J^s(0) = [\xi_1\ \xi_2\ \cdots\ \xi_n]$.** — **True**

Every exponential is I at θ = 0, and $\mathrm{Ad}_I$ is the identity, so each column is just the reference twist. This is why "assume θ = 0" makes velocity and torque problems quick.

*Body version:* $J^b(0) = \mathrm{Ad}^{-1}_{g_{st}(0)}[\xi_1\cdots\xi_n]$: the same axes, re‑measured from the tool origin.

---

### 14.6 The patterns behind all 47

| Hidden condition | Statements it decides |
|---|---|
| $\Vert\omega\Vert= 1$ required | Rodrigues (6), unit twist (19), screw exp (20) |
| ω ≠ 0 vs. ω = 0 | exp formula (20), pitch ∞, prismatic twist |
| Pitch 0 vs. ≠ 0 | $v = -\omega\times q$ (17), fixed axis points (34) |
| "for all" vs. "there exists" | 6, 9, 35 |
| Which frame or coordinates | 1, 28, 30, 39, 40, 46 |
| Which group operation | 2, 3 |
| Surjective vs. injective | 7, 8, 13, 22 |
| Commuting vs. not | 4, 24, 25, 37 |
| Point vs. vector (bottom 1 vs. 0) | 15, 16, 21, 29 |
| Jacobian dependence | 43, 44, 45, 47 |

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
