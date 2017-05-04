import matplotlib.pyplot as plt
import numpy as np
import tensorflow as tf

from util import DataReader, CrossValidation, partition_data


fp = '../data/E-GEOD-48350/E-GEOD-48350-combined.csv'


def FFNN(x, weights, biases, keep_prob=0.5, h_l=2):
    h = tf.add(tf.matmul(x, weights['hidden_1']), biases['hidden_1'])
    h = tf.nn.dropout(tf.tanh(h), 0.5)
    if h_l == 2:
        h = tf.add(tf.matmul(h, weights['hidden_2']), biases['hidden_2'])
        h = tf.nn.dropout(tf.tanh(h), 0.5)
    out = tf.add(tf.matmul(h, weights['out']), biases['out'])
    return out

def MAX(x, weights, biases, dim):
    return tf.concat([tf.ones([dim, 1], dtype=tf.float64), tf.zeros([dim, 1], dtype=tf.float64)], 1)

# Parameters
n_hidden_1 = 2048
n_hidden_2 = 128
n_classes = 2
fold = 10 # Cross validation
learning_rate = 0.01
dropout_prob = 0.5 # keep_prob = 1 - dropout_prob

# Read and parse input data
x, y = DataReader(fp).get_data()
partition = partition_data(x, y, [0.8, 0.2])
data = CrossValidation(partition[0][0], partition[0][1], fold)
test = np.array(partition[1][0]), np.array(partition[1][1])
n_input = len(x[0])

# Tensorflow variables
weights = {
    'hidden_1': tf.Variable(tf.random_normal([n_input, n_hidden_1])),
    'hidden_2': tf.Variable(tf.random_normal([n_hidden_1, n_hidden_2])),
    'out': tf.Variable(tf.random_normal([n_hidden_2, n_classes]))
}

biases = {
    'hidden_1': tf.Variable(tf.random_normal([n_hidden_1])),
    'hidden_2': tf.Variable(tf.random_normal([n_hidden_2])),
    'out': tf.Variable(tf.random_normal([n_classes]))
}

x = tf.placeholder(tf.float32, [None, n_input])
y = tf.placeholder(tf.float32, [None, n_classes])
dim = tf.placeholder(tf.int32)

model = FFNN #either MAX or FFNN

if model == MAX:
    pred = model(x, weights, biases, dim)
elif model == FFNN:
    pred = model(x, weights, biases, keep_prob=1-dropout_prob)
cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=pred, labels=y))
if model != MAX:
    optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(cost)
correct_pred = tf.equal(tf.argmax(pred, 1), tf.argmax(y, 1))
accuracy = tf.reduce_mean(tf.cast(correct_pred, tf.float32))

init = tf.global_variables_initializer()

v_lss_hist = []
t_lss_hist = []

# Train
with tf.Session() as sess:
    sess.run(init)
    steps = 2001
    print('Iteration | Valid loss | Valid acc  | Test loss  | Test acc   ')
    for i in range(steps):
        train_x, train_y = data.get_train_set()
        sess.run(optimizer, feed_dict={x: train_x, y: train_y})

        test_x, test_y = data.get_train_set()
        v_acc = sess.run(accuracy, feed_dict={x: test_x, y: test_y, dim: len(test_x)})
        v_lss = sess.run(cost, feed_dict={x: test_x, y: test_y, dim: len(test_x)})
        t_acc = sess.run(accuracy, feed_dict={x: test[0], y: test[1], dim: len(test[0])})
        t_lss = sess.run(cost, feed_dict={x: test[0], y: test[1], dim: len(test[0])})
        v_lss_hist.append(v_lss)
        t_lss_hist.append(t_lss)
        print('{0:9} |   {1:.6f} |   {2:.6f} |   {3:.6f} |    {4:.6f}'.format(i, v_lss, v_acc, t_lss, t_acc))


# Plot validation set and test set loss over epochs
x_r = range(0, 2001, 10)
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.scatter(x_r, v_lss_hist, label='Validation loss')
ax1.scatter(x_r, t_lss_hist, label='Test loss')
plt.legend(loc='upper left')
plt.show()