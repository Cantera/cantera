# Writing Tests

## Python Tests

- To allow use of deprecated methods in Python tests, use the decorator:
  ```python
  @pytest.mark.usefixtures("allow_deprecated")
  ```

### Transition from `unittest` to `pytest`

Starting in Cantera 2.6, we transitioned from using the built-in `unittest` module to
run Cantera's Python test suite to using the `pytest` framework. However, many of the
existing tests are still written in the `unittest` style (which is also supported by
`pytest`). The following recommendations apply to new tests and tests that are being
significantly modified:

- Use `pytest`-style asserts instead of `unittest` member functions like
  `assertWhatever`. Examples:
  ```python
  assert 0 <= idom < len(self.domains)
  assert gas.n_species == N + 1
  assert "X_H2" in header.split(",")
  ```

- Use `pytest.approx` in preference to `assertAlmostEqual`. Add
  `from pytest import approx` to any test file that doesn't already contain it.
  Examples:
  ```python
  assert q1.T == approx(T1)  # scalar
  assert sl_mole == approx(sl_mass, rel=0.1)  # specific relative tolerance
  assert list(jet.X[k, :]) == approx(list(self.sim.X[k, :]), tol_X)  # arrays
  assert X == approx([0.07748481, 0.048165072, 0.01446654])
  ```

- Use `with pytest.raises` instead of `unittest`'s method of handling "assert raises".
  Example:
  ```python
  with pytest.raises(ct.CanteraError, match="No component named 'spam'"):
      self.r1.set_advance_limit("spam", 0.1)
  ```
