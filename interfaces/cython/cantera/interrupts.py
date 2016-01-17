def no_op(t):
    """
    This function does nothing. It is used as an interrupt in the 1D solver
    C++ loop where a pure Python function is needed in order for
    KeyboardInterrupt events to be captured.
    """
    return 0.0
