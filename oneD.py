from typing import Union

import numpy as np
from typeguard import typechecked

from constants import ureg


class Operation(object):
    def __init__(self, *operands):
        self.operands = operands

    def compute(self, method=None):
        raise NotImplementedError


class UnaryOperation(Operation):
    def __init__(self, a):
        super().__init__(a)
        self._a = a

    @property
    def a(self):
        return self._a

    @a.setter
    def a(self, a):
        self._a = a

    def compute(self, method=None):
        raise NotImplementedError


class BinaryOperation(Operation):
    def __init__(self, a, b):
        super().__init__(a, b)
        self._a = a
        self._b = b

    @property
    def a(self):
        return self._a

    @a.setter
    def a(self, a):
        self._a = a

    @property
    def b(self):
        return self._b

    @b.setter
    def b(self, b):
        self._b = b

    def compute(self, method=None):
        raise NotImplementedError


class Value(UnaryOperation):
    def __init__(self, value):
        super().__init__(value)
        self.value = self.a

    def compute(self, method=None):
        return self

    def __truediv__(self, other):
        return Value(self.a / other)

    def __rtruediv__(self, other):
        return Value(other / self.a)


class Variable(Value):
    def __init__(self, value):
        super().__init__(value)
        self.value = self.a


class Differential(Value):
    def __init__(self, value: Variable):
        super().__init__(value)
        self.value = self.a


class Partial(BinaryOperation):
    def __init__(self, of, respect: Differential):
        super().__init__(of, respect)
        self.of = self.a
        self.respect = self.b

    def compute(self, method=None):
        raise NotImplementedError


class TensorField(Value):
    def __init__(self, value):
        super().__init__(value)


class VectorField(TensorField):
    def __init__(self, value):
        super().__init__(value)


class ScalarField(TensorField):
    def __init__(self, value):
        super().__init__(value)


class Scalar(Value):
    def __init__(self, value):
        super().__init__(value)


class Constant(Value):
    def __init__(self, value):
        super().__init__(value)


class Gradient(UnaryOperation):
    @typechecked
    def __init__(self, scalar: ScalarField):
        super().__init__(scalar)
        self.scalar = self.a

    def compute(self, method=None) -> VectorField:
        raise NotImplementedError


class Divergence(UnaryOperation):
    @typechecked
    def __init__(self, vector: VectorField):
        super().__init__(vector)
        self.vector = self.a

    def compute(self, method=None) -> ScalarField:
        raise NotImplementedError


class MaterialDerivative(BinaryOperation):
    @typechecked
    def __init__(self, tensor: TensorField, differential: Differential):
        super().__init__(tensor, differential)
        self.tensor = self.a
        self.scalar = self.b

    def compute(self, method=None):
        raise NotImplementedError


class Equation(object):
    def __init__(self, left, right):
        self.left = left
        self.right = right



