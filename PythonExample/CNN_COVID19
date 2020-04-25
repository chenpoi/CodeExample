import pandas as pd
from PIL import Image
import numpy as np
import tensorflow as tf
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
from collections import Counter


def plot_confusion_matrix(cm, classes,
                          normalize=False,
                          title='Confusion matrix',
                          cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    import itertools
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    print(cm)

    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45)
    plt.yticks(tick_marks, classes)

    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, format(cm[i, j], fmt),
                 horizontalalignment="center",
                 color="white" if cm[i, j] > thresh else "black")

    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    plt.tight_layout()


meta_df = pd.read_csv("./Chest_xray_Corona_Metadata.csv")
meta_df["path"] = \
    ["./Coronahack-Chest-Xray-Dataset/Coronahack-Chest-Xray-Dataset/" +
     meta_df.loc[i, "Dataset_type"].lower() +
     "/" +
     meta_df.loc[i, "X_ray_image_name"]
     for i in meta_df.index]

meta_df["resized_path"] = \
    ["./Coronahack-Chest-Xray-Dataset/resized_image/" +
     meta_df.loc[i, "X_ray_image_name"]
     for i in meta_df.index]

meta_df["final_lable"] = \
    [meta_df.loc[i, "Label"] +
     "_" +
     str(meta_df.loc[i, "Label_2_Virus_category"]) +
     "_" +
     str(meta_df.loc[i, "Label_1_Virus_category"])
     for i in meta_df.index]

print(meta_df["final_lable"].unique())
meta_df = meta_df.loc[[i in ["Normal_nan_nan", "Pnemonia_nan_bacteria", "Pnemonia_COVID-19_Virus", "Pnemonia_SARS_Virus"]
                       for i in meta_df["final_lable"]], ]


all_label = {}
for cnt, i in enumerate(meta_df["final_lable"].unique()):
    temp = np.zeros(len(meta_df["final_lable"].unique()))
    temp[cnt] = 1
    print(i, temp)
    all_label[i] = temp

n_class = len(all_label)


# 1. resize the image

resized_H = 150
resized_W = 150

print(len(meta_df))
for i in meta_df.index:
    my_image = Image.open(meta_df.loc[i, "path"]).convert('L')
    new_image = my_image.resize((resized_H, resized_W))
    new_image.save(meta_df.loc[i, "resized_path"])

# 1.5 small clsuter oversampling
label_counter = Counter(meta_df["final_lable"])
print(label_counter)
sampling_target = label_counter.most_common()[0][1]
for label in label_counter.keys():
    print(label)
    temp_df = meta_df.loc[meta_df["final_lable"] == label, ]

    multiply_n = min(sampling_target // label_counter[label] - 1, 50)
    addition_n = sampling_target % label_counter[label]
    final_df = pd.DataFrame()
    if multiply_n:
        final_df = pd.concat([temp_df] * multiply_n, ignore_index=True)
    if addition_n:
        final_df = pd.concat([final_df, temp_df.iloc[:addition_n,]])
    if multiply_n or addition_n:
        meta_df = pd.concat([meta_df, final_df])

print(meta_df.shape)

# 2. train, test split. (X and Y)

Xtrain, Xtest, Ytrain, Ytest = train_test_split(
    meta_df.loc[:, "resized_path"],
    np.array([all_label[i] for i in meta_df.loc[:, "final_lable"]]).reshape([-1, n_class]), test_size=0.33)

image_test = np.array([np.array(Image.open(Xtest.iloc[i])).reshape(
    [resized_H * resized_W]
)
    for i in range(len(Xtest))])



# 3. prepare convolution layer functions


def init_weights(shape):
    init_random_dist = tf.truncated_normal(shape, stddev=0.1)
    return tf.Variable(init_random_dist)


def init_bias(shape):
    init_bias_vals = tf.constant(0.1, shape=shape)
    return tf.Variable(init_bias_vals)


def conv2d(x, W):
    return tf.nn.conv2d(x, W, strides=[1, 1, 1, 1], padding='SAME')


def max_pool_2by2(x):
    return tf.nn.max_pool(x, ksize=[1, 2, 2, 1],
                          strides=[1, 2, 2, 1], padding='SAME')


def convolutional_layer(input_x, shape):
    W = init_weights(shape)
    b = init_bias([shape[3]])
    return tf.nn.relu(conv2d(input_x, W) + b)


def normal_full_layer(input_layer, size):
    input_size = int(input_layer.get_shape()[1])
    W = init_weights([input_size, size])
    b = init_bias([size])
    return tf.matmul(input_layer, W) + b


def duplicated_continuous_larers(input_layer, shape, n_duplicated):
    if n_duplicated:
        convo = convolutional_layer(input_layer, shape=shape)
        shape[2] = shape[3]
        return duplicated_continuous_larers(convo, shape, n_duplicated - 1)
    else:
        return input_layer


# 4. build graph and declare tensorflow variables


x = tf.placeholder(tf.float32, shape=[None, resized_H * resized_W])
y_true = tf.placeholder(tf.float32, shape=[None, n_class])
x_image = tf.reshape(x, [-1, resized_H, resized_W, 1])

convo_1 = duplicated_continuous_larers(x_image, shape=[3, 3, 1, 32], n_duplicated=1)
convo_1_pooling = max_pool_2by2(convo_1)

convo_2 = duplicated_continuous_larers(convo_1_pooling, shape=[3, 3, 32, 64], n_duplicated=1)
convo_2_pooling = max_pool_2by2(convo_2)

convo_3 = duplicated_continuous_larers(convo_2_pooling, shape=[3, 3, 64, 64], n_duplicated=1)
convo_3_pooling = max_pool_2by2(convo_3)

convo_2_flat = tf.reshape(convo_3_pooling, [-1, 23104])

full_layer_one = tf.nn.relu(normal_full_layer(convo_2_flat, 512))

hold_prob = tf.placeholder(tf.float32)
full_one_dropout = tf.nn.dropout(full_layer_one, keep_prob=hold_prob)

y_pred = normal_full_layer(full_one_dropout, n_class)

cross_entropy = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(labels=y_true, logits=y_pred))
optimizer = tf.train.AdamOptimizer(learning_rate=0.0001)
train = optimizer.minimize(cross_entropy)


# 5. run session


def next_batch(num, data, labels):
    # Return a total of `num` random samples and labels.
    data.index = range(len(data))
    idx = np.arange(0, len(data))
    np.random.shuffle(idx)
    idx = idx[:num]
    data_shuffle = [np.array(Image.open(data[i])).reshape([resized_H * resized_W]) for i in idx]
    labels_shuffle = [labels[i] for i in idx]

    return np.asarray(data_shuffle), np.asarray(labels_shuffle).reshape(num, n_class)


init = tf.global_variables_initializer()
steps = 500
figures = []
with tf.Session() as sess:
    sess.run(init)

    for i in range(steps):

        batch_x, batch_y = next_batch(49, Xtrain, Ytrain)

        sess.run(train, feed_dict={x: batch_x, y_true: batch_y, hold_prob: 0.5})

        # PRINT OUT A MESSAGE EVERY 100 STEPS
        if i % 100 == 0:
            print('Currently on step {}'.format(i))
            print('Accuracy is:')
            # Test the Train Model
            matches = tf.equal(tf.argmax(y_pred, 1), tf.argmax(y_true, 1))
            acc = tf.reduce_mean(tf.cast(matches, tf.float32))

            cnf_matrix = tf.confusion_matrix(
                labels=tf.argmax(y_pred, 1),
                predictions=tf.argmax(y_true, 1),
                num_classes=n_class
            )
            cnf_res = sess.run(cnf_matrix, feed_dict={x: image_test, y_true: Ytest, hold_prob: 1.0})
            fig = plt.figure()
            plot_confusion_matrix(cnf_res, classes=all_label.keys(),
                                  title='Confusion matrix, without normalization')
            figures.append(fig)

            print(sess.run(acc, feed_dict={x: image_test, y_true: Ytest, hold_prob: 1.0}))

            print('\n')

# 6. print result and quality
fig = plt.figure()
figures[0]
figures[1]
figures[2]
figures[3]
plt.show()
