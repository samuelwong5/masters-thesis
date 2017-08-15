import matplotlib.pyplot as plt
import numpy as np
import tensorflow as tf

from util import DataReader, CrossValidation, partition_data


#fp = '../data/E-GEOD-48350/E-GEOD-48350-combined.csv'
fp = '../data/48350-AD_6-HD.csv'

def FFNN(x, weights, biases, keep_prob=0.5, layers=1):
    h = x
    for i in range(layers):
        layer_name = 'hidden_{}'.format(i+1)
        h = tf.add(tf.matmul(h, weights[layer_name]), biases[layer_name])
        h = tf.nn.dropout(tf.tanh(h), keep_prob)
    out = tf.add(tf.matmul(h, weights['out']), biases['out'])
    return out

def AutoEncoder(x, weights, biases, act=tf.nn.sigmoid):
    h = act(tf.add(tf.matmul(x, weights['encode_1']), biases['encode_1']))
    h = act(tf.add(tf.matmul(h, weights['encode_2']), biases['encode_2']))
    h = act(tf.add(tf.matmul(h, weights['decode_1']), biases['decode_1']))
    h = act(tf.add(tf.matmul(h, weights['decode_2']), biases['decode_2']))
    return h

# Autoencoder params
autoencoder_1 = 1024
autoencoder_2 = 248

# Parameters
n_hidden_3 = 128
n_classes = 2
fold = 10 # Cross validation
learning_rate = 1
dropout_prob = 0.5 # keep_prob = 1 - dropout_prob

# Read and parse input data
x, y = DataReader(fp).get_data()

def inverse_argmax(l):
    m = len(set(l))
    for i in range(len(l)):
        hold = l[i]
        l[i] = [0] * m
        l[i][hold] = 1

inverse_argmax(y)

partition = partition_data(x, y, [0.8, 0.2])
data = CrossValidation(partition[0][0], partition[0][1], fold)
test = np.array(partition[1][0]), np.array(partition[1][1])
n_input = len(x[0])

autoencoder_encode_weights_1 = tf.Variable(tf.random_normal([n_input, autoencoder_1]))
autoencoder_encode_weights_2 = tf.Variable(tf.random_normal([autoencoder_1, autoencoder_2]))
autoencoder_encode_biases_1 = tf.Variable(tf.random_normal([autoencoder_1]))
autoencoder_encode_biases_2 = tf.Variable(tf.random_normal([autoencoder_2]))


ae_biases = {
    'encode_1': autoencoder_encode_biases_1,
    'encode_2': autoencoder_encode_biases_2,
    'decode_1': tf.Variable(tf.random_normal([autoencoder_1])),
    'decode_2': tf.Variable(tf.random_normal([n_input]))
}

ae_weights = {
    'encode_1': autoencoder_encode_weights_1,
    'encode_2': autoencoder_encode_weights_2,
    'decode_1': tf.Variable(tf.random_normal([autoencoder_2, autoencoder_1])),
    'decode_2': tf.Variable(tf.random_normal([autoencoder_1, n_input]))
}

# Tensorflow variables
weights = {
    'hidden_1': autoencoder_encode_weights_1,
    'hidden_2': autoencoder_encode_weights_2,
    'hidden_3': tf.Variable(tf.random_normal([autoencoder_2, n_hidden_3])),
    'out': tf.Variable(tf.random_normal([n_hidden_3, n_classes]))
}

biases = {
    'hidden_1': autoencoder_encode_biases_1,
    'hidden_2': autoencoder_encode_biases_2,
    'hidden_3': tf.Variable(tf.random_normal([n_hidden_3])),
    'out': tf.Variable(tf.random_normal([n_classes]))
}

x = tf.placeholder(tf.float32, [None, n_input])
y = tf.placeholder(tf.float32, [None, n_classes])
dim = tf.placeholder(tf.int32)


autoencoder_model = AutoEncoder(x, ae_weights, ae_biases)
autoencoder_cost = tf.reduce_mean(tf.pow(x - autoencoder_model, 2))
autoencoder_optimizer = tf.train.AdamOptimizer(learning_rate=0.5).minimize(autoencoder_cost)

pred = FFNN(x, weights, biases, keep_prob=1-dropout_prob, layers=3)
cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=pred, labels=y))
global_step = tf.Variable(0, trainable=False)
starter_learning_rate = 1.0
learning_rate = tf.train.exponential_decay(starter_learning_rate, global_step,
                                           50, 0.5, staircase=False)
optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(cost)

correct_pred = tf.equal(tf.argmax(pred, 1), tf.argmax(y, 1))
accuracy = tf.reduce_mean(tf.cast(correct_pred, tf.float32))

init = tf.global_variables_initializer()

t_acc_hist = []
v_lss_hist = []
t_lss_hist = []

# Train
with tf.Session() as sess:
    sess.run(init)
    print('Training autoencoder...')
    steps = 200
    for i in range(steps):
        train_x, _ = data.get_train_set()
        sess.run(autoencoder_optimizer, feed_dict={x: train_x})
        test_x, _ = data.get_train_set()
        v_lss = sess.run(autoencoder_cost, feed_dict={x: test_x})
        print('{0:9} | {1:2.5f}'.format(i+1, v_lss))

    print('Training classifier..')
    steps = 200
    print('Iteration | Valid loss | Valid acc  | Test loss  | Test acc   ')
    for i in range(steps):
        train_x, train_y = data.get_train_set()
        sess.run(optimizer, feed_dict={x: train_x, y: train_y, global_step: i})

        test_x, test_y = data.get_train_set()
        v_acc = sess.run(accuracy, feed_dict={x: test_x, y: test_y, dim: len(test_x)})
        v_lss = sess.run(cost, feed_dict={x: test_x, y: test_y, dim: len(test_x)})
        t_acc = sess.run(accuracy, feed_dict={x: test[0], y: test[1], dim: len(test[0])})
        t_lss = sess.run(cost, feed_dict={x: test[0], y: test[1], dim: len(test[0])})
        #v_lss_hist.append(v_lss)
        t_acc_hist.append(t_acc)
        #t_lss_hist.append(t_lss)
        #print('{0:9} |   {1:2.5f} |   {2:2.5f} |   {3:2.5f} |    {4:2.4f}'.format(i, v_lss, v_acc, t_lss, t_acc))



# Plot validation set and test set loss over epochs
'''
x_r = range(0, 200)
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.scatter(x_r, v_lss_hist, label='Validation loss')
ax1.scatter(x_r, t_lss_hist, label='Test loss')
plt.legend(loc='upper left')
plt.show()
'''
print(sum(t_acc_hist[-10:]) / 10)
print(len(test[0]))