import numpy as np
from sklearn.datasets import load_breast_cancer

from sklearn import svm


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

class my_scaler:

    def __int__(self, X=None):
        self.data = X

    def fit(self, X):
        self.mean = np.mean(X, axis=0)
        self.std = np.std(X, axis=0)

    def transform(self, X):
        if np.ndim(X) - np.ndim(self.mean) == 1:
            return (X - self.mean)/self.std
        else:
            ValueError("the dimension of train data don't match the transformed data")


def split_train_test(X, Y, test_size=0.5, psudo_seed=1):
    test_n = round(test_size * X.shape[0])

    if X.shape[0] == Y.shape[0]:
        np.random.seed(psudo_seed)
        np.random.shuffle(X)
        np.random.seed(psudo_seed)
        np.random.shuffle(Y)

        test_X = X[:test_n]
        test_Y = Y[:test_n]
        train_X = X[test_n:]
        train_Y = Y[test_n:]

        return train_X, test_X, train_Y, test_Y
    else:
        ValueError("X and Y should have same number of records")




def medical():
  data = load_breast_cancer()
  X, Y = data.data, data.target
  Xtrain, Xtest, Ytrain, Ytest = split_train_test(X, Y, test_size=0.33)
  return Xtrain, Xtest, Ytrain, Ytest, rbf, 1e-3, 200


Xtrain, Xtest, Ytrain, Ytest, kernel, lr, n_iters = medical()
# make sure the targets are (-1, +1)
Ytrain[Ytrain == 0] = -1
Ytest[Ytest == 0] = -1

scaler = my_scaler()
scaler.fit(Xtrain)

Xtrain = scaler.transform(Xtrain)
Xtest = scaler.transform(Xtest)


 # now we'll use our custom implementation
model = mySvm(kernel=kernel, C=1.0)

model.fit(Xtrain, Ytrain, learn_rate=lr, n_perm=n_iters)
print("precisions is " + str(model.precisions(Xtest, Ytest)))


###################################################
# compare it to SVC in sklearn
###################################################


sk_svc = svm.SVC(kernel=kernel, C=1.0)
sk_svc.fit(Xtrain, Ytrain)
print("precisions in sklearn is " + str(sk_svc.score(Xtest, Ytest)))
