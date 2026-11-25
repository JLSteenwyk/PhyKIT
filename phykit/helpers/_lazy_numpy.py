"""Defer NumPy imports while keeping attribute caches local to each caller."""


class LazyNumpy:
    _module = None

    def __getattr__(self, name):
        module = self._module
        if module is None:
            import numpy

            module = self._module = numpy

        value = getattr(module, name)
        setattr(self, name, value)
        return value
