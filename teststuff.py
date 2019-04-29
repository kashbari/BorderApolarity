def partnk(n, k):
    def _partnk(n, k, pre):
        if n <= 0:
            return []
        if k == 1:
            if n <= pre:
                return [[n]]
            return []
        ret = []
        for i in range(min(pre, n), 0, -1):
            ret += [[i] + sub for sub in _partnk(n-i, k-1, i)]
        return ret
    return _partnk(n, k, n)
