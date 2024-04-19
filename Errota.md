**Updates on 19 April 2024**

In my codes, there exists an error in the statement:

```
X = X + Sigma_x*H'*S_square*phi(S_square*(Y(i) - H*X), k);
```

Please change to:

```
X = X + Sigma_x*H'*S_square'*phi(S_square*(Y(i) - H*X), k);
```

That is, applying a transpose symbol on the first variable "S_square".

This is because

```
S^(-1) = S_square' * S_square;    % By Cholesky Factorization: S_square = chol(S^(-1))
```

This typo did not introduce fatal errors to the experiments in my paper because the variable "S_square" therein is one-dimensional. However, this will cause potential errors for multi-dimensional cases.
