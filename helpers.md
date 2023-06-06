# Some helper functions for Liam Eagen's protocol.

Here, we assume that we are over a curve $E$, defined over a field $F$ in Weierstrass form (and let's also assume the equation has the form $y^2 - x^3 - c = 0$)

Let's denote the equation $\text{EQ}=0$.


## Normal form

For a polynomial of two variables $P(x,y)$ let's call its *normal form* a canonical representative of $P(x,y) \mod \text{EQ}$, obtained by replacing $y^2 \rightarrow x^3 + c$, until it is no longer possible.

### Algorithm:

Represent polynomial as $p_{i}(x)y^i$. Compute the sum:
$\text{nf}(P)=p_i(x)(x^3+c)^{i \text{//} 2}y^{i\%2}$.

**Definition:** regular functions are polynomials (mod EQ).

We can implement multiplication of regular functions by first multiplying them as polynomials, and then computing the normal form.

---

Our main goal is constructing a regular function which vanishes exactly in the set of points $A_1, .., A_n \in E$. We assume that these points are nonzero (and if any of them are zero, they should be universally ignored).

Such function only exists iff sum of these points is 0 (in a sense of elliptic curve addition).

Such function has degree $n$, which means that it has highest coefficient of degree $n$, with $x$ treated as having degree $2$ and $y$ treated as having degree $3$.

(for example:

n=2 - function has the form $a+bx$

n=3 - function has the form $a+bx+cy$

n=4 - function has the form $a+bx+cx^2 + dy$

et cetera)

**Definition:** point-set witness is the following structure:

```
inputs: Vec<E>
output: E
f: RegularFunction
```

Where ``f`` vanishes only in elements of input and in the output (=> ``output == -sum(inputs)``).

The rationale for separating them is the idea that such structures can be "merged" together.

We define the following binary "merge" method:

```
fn merge(a: &PointSetWitness, b: &PointSetWitness) -> PointSetWitness
```

On inputs and outputs it should work in a straightforward way required by:

```
merge(a,b).inputs == a.inputs.clone().append(b.inputs());
merge(a,b).output == a.outputs + b.output
```

Now, let us explain how we obtain ``merge.f``. First of all, if either of ``a.output`` or ``b.output`` are 0, we can just set ``merge(a,b).f == a.f * b.f``.


In what follows, we will use ``linefunc``, satisfying following requirements:
<details><summary>linefunc requirements:</summary>

```
linefunc(a: E, b: E) -> RegularFunction
```

1. Must throw if both $a$ and $b$ are $0$.
2. Must correctly return equation of a line passing through a pair of points if both are nonzero and unequal.
3. Must correctly return the equation of a tangent line if points coincide.
4. What it does if one of the points is zero, and other is not is left up to the implementer (could return vertical line, or throw).
</details>

<br>
Now, let's compute:

```
let numerator = a.f * b.f * linefunc(-a.output, -b.output)
```

vanishes in all the inputs, required output (``a.output + b.output``), and following four points that we do not need: ``a.output``, ``b.output``, ``-a.output``, ``-b.output``. Remember - computations are done modulo ``EQ`` here and everywhere later.

The plan is to divide numerator by two vertical lines: ``linefunc(a.output, -a.output)`` and ``linefunc(b.output, -b.output)``.

For this, we need a partial evaluation function (which I will call "subst" and hope its syntax can be inferred).

```
let num_restricted = numerator.subst(Var::x, a.output.x)
let eq_restricted = EQ.subst(Var::x, a.output.x)
let quot = num_restricted / eq_restricted // this is division of univariate polynomials, moreover, eq_restricted has the form y^2 + const

let numerator = (numerator - RegularFunction(quot) * EQ)/linefunc(a.output, -a.output) // the division by the vertical line is essentially a division of univariate polynomial (by x-const)

// what happened from math standpoint is that we have ensured that (numerator - smth*EQ) vanishes on this vertical line, and now can divide

[...] repeat the same once again for b.output
```

And this is how we compute ``merge(a, b).f``.

## Computing the witness for a collection of points.

Assume we are given the collection of points ``A[0], ..., A[n-1]``. Initialize the ``Q[i]: PointSetWitness`` using the following initial data:

```
Q[i].inputs = vec![A[i]]
Q[i].output = -A[i]
Q[i].f = linefunc(A[i], -A[i])
```

And now merge them (organized in a binary tree). The final merge's ``output`` must be zero (or function throws because the initial sum, evidently, was nonzero), and the function should return final merge's ``f``.

There are few minor optimizations - for example, if we have obtained sum subtree which already got output == 0, we do not need to merge it, we can just push the corresponding polynomial f to some temporary array, and then multiply all the polynomials in this array.

We should just drop all the zero points before we even started the process.

The function is subject to the following tests:

1. The output function f=a(x)+y*b(x) vanishes in all the points ``A[0], ..., A[n-1]``
   
2. a(x) has degree k//2, and b(x) has degree (k-1)//2, where k is amount of nonzero points in A[0], ..., A[n-1]. In particular, that means that the total potential amount of coefficients is k//2 + (k-1)//2 + 2 = k.
   
3. Function works correctly for arrays with multiple zeros, and multiple negations / low scalar multiples of the same point.

4. Norm of a function (defined as ``f(x,y)*f(x,-y)==a**2 - (x**3 + c)*b**2``, normalized to have highest coefficient == 1, equals product of ``x-A[i].x``).