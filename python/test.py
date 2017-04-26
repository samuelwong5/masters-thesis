import tensorflow as tf
import random
import numpy as np

learning_rate = 0.05
batch_size = 128

n_input = 1
n_steps = 50
n_hidden = 64
n_classes = 2

x = tf.placeholder(tf.float32, [None, n_steps, n_input])
y = tf.placeholder(tf.float32, [None, n_classes])

weights = {
    'out': tf.Variable(tf.random_normal([n_hidden, n_classes]))
}

biases = {
    'out': tf.Variable(tf.random_normal([n_classes]))
}

def RNN(x, weights, biases):
    #x = tf.unstack(x, n_steps, 1)

    lstm_cell = tf.contrib.rnn.LSTMCell(n_hidden)

    outputs, states = tf.nn.dynamic_rnn(lstm_cell, x, dtype=tf.float32)

    outputs = tf.stack(outputs)
    outputs = tf.transpose(outputs, [1, 0, 2])

    return tf.matmul(outputs[-1], weights['out']) + biases['out']


def SequenceBatchGen(batch_size):
    x = []
    y = []
    for i in range(batch_size):
        if random.randint(0, 1) == 0:
            x.append(list(range(n_steps)))
            y.append([1,0])
        else:
            x.append([random.randint(0, 50) for i in range(n_steps)])
            y.append([0,1])
    return x, y


pred = RNN(x, weights, biases)

cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=pred, labels=y))

optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(cost)

correct_pred = tf.equal(tf.argmax(pred, 1), tf.argmax(y, 1))

accuracy = tf.reduce_mean(tf.cast(correct_pred, tf.float32))

init = tf.global_variables_initializer()

with tf.Session() as sess:
    sess.run(init)
    steps = 10001
    for i in range(steps):
        b_x, b_y = SequenceBatchGen(batch_size)
        b_x = np.array(b_x)
        b_y = np.array(b_y)
        b_x = b_x.reshape((batch_size, n_steps, n_input))
        sess.run(optimizer, feed_dict={x: b_x, y: b_y})
        if i % 100 == 0:
            acc = sess.run(accuracy, feed_dict={x: b_x, y: b_y})

            loss = sess.run(cost, feed_dict={x: b_x, y: b_y})

            print('Iteration {0}, Loss: {1:.6f}, Acc: {2:.5f}'.format(i, loss, acc))

