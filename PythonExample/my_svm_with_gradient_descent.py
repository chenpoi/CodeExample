import numpy as np
from sklearn.datasets import load_breast_cancer
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler


def linear(X1, X2, c=0):
    return X1.dot(X2.T) + c


def rbf(X1, X2, gamma=None):
    if gamma is None:
        gamma = 1.0 / X1.shape[-1] # 1 / D
      # gamma = 0.05
      # gamma = 5. # for donut and spiral
    if np.ndim(X1) == 1 and np.ndim(X2) == 1:
        result = np.exp(-gamma * np.linalg.norm(X1 - X2)**2)
    elif (np.ndim(X1) > 1 and np.ndim(X2) == 1) or (np.ndim(X1) == 1 and np.ndim(X2) > 1):
        result = np.exp(-gamma * np.linalg.norm(X1 - X2, axis=1)**2)
    elif np.ndim(X1) > 1 and np.ndim(X2) > 1:
        result = np.exp(-gamma * np.linalg.norm(X1[:, np.newaxis] - X2[np.newaxis, :], axis=2)**2)
    return result


def sigmoid(X1, X2, gamma=0.05, c=1):
    return np.tanh(gamma * X1.dot(X2.T) + c)


class mySvm:

    def __init__(self, C, kernel):
        self.C = C
        self.kernel = kernel


    def _objective_function(self, XYY):
        return np.sum(self.alpha) - 0.5 * np.sum(np.outer(self.alpha, self.alpha) * XYY)


    def fit(self, X, Y, learn_rate=1e-5, n_perm=400):
        n_samples, n_features = X.shape
        losses = []

        self.alpha = np.zeros(n_samples)
        self.X = X
        self.Y = Y


        K = self.kernel(X, X)
        YY = np.outer(Y, Y)
        XYY = K * (YY)


        for _ in range(n_perm):
            loss = self._objective_function(XYY)
            losses.append(loss)

            gr = np.ones(n_samples) - self.alpha.dot(XYY)
            self.alpha = self.alpha + learn_rate * gr

            self.alpha[self.alpha < 0] = 0
            self.alpha[self.alpha > self.C] = self.C

        self.bs = np.mean(Y[self.alpha == 0] - ((self.alpha * Y).dot(self.kernel(X, X[self.alpha == 0]))))


    def _decision_function(self, X):
        return (self.alpha * self.Y).dot(self.kernel(self.X, X)) + self.bs

    def predict(self, X):
        return np.sign(self._decision_function(X))

    def precisions(self, X, Y):
        P = self.predict(X)
        return np.mean(Y == P)


def medical():
  data = load_breast_cancer()
  X, Y = data.data, data.target
  Xtrain, Xtest, Ytrain, Ytest = train_test_split(X, Y, test_size=0.33)
  return Xtrain, Xtest, Ytrain, Ytest, rbf, 1e-3, 200


Xtrain, Xtest, Ytrain, Ytest, kernel, lr, n_iters = medical()
# make sure the targets are (-1, +1)
Ytrain[Ytrain == 0] = -1
Ytest[Ytest == 0] = -1

  # scale the data
scaler = StandardScaler()
Xtrain = scaler.fit_transform(Xtrain)
Xtest = scaler.transform(Xtest)

  # now we'll use our custom implementation
model = mySvm(kernel=kernel, C=1.0)

model.fit(Xtrain, Ytrain, learn_rate=lr, n_perm=n_iters)
print("preciisions is " + str(model.precisions(Xtest, Ytest)))
