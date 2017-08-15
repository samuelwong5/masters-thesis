import matplotlib.pyplot as plt
import numpy as np
import tensorflow as tf

from util import DataReader, CrossValidation, partition_data


# Each task has individual layer which is concatenated with shared layer.
def mtl_2(x, weights, biases, keep_prob=0.5, act=tf.tanh):
    h_shared = tf.add(tf.matmul(x, weights['hidden_shared']), biases['hidden_shared'])
    h_task = tf.add(tf.matmul(x, weights['hidden_task']), biases['hidden_task'])
    h_shared = tf.nn.dropout(act(h_shared), keep_prob)
    h_task = tf.nn.dropout(act(h_task), keep_prob)


    h_concat = tf.concat([h_shared, h_task], 1)
    h = tf.add(tf.matmul(h_concat, weights['hidden_2']), biases['hidden_2'])
    h = tf.nn.dropout(act(h), keep_prob)


    out = tf.add(tf.matmul(h, weights['out']), biases['out'])
    return out


def gridsearch(params):
    # Parameters
    n_hidden_shared = params['n_hidden_shared']
    n_hidden_task = params['n_hidden_task']
    n_hidden_2 = params['n_hidden_2']
    act_func = params['act_func']
    n_classes = 2
    fold = 10 # Cross validation
    learning_rate = 1
    dropout_prob = 0.5 # keep_prob = 1 - dropout_prob

    ad_data = '../data/48350-AD_6-HD.csv'
    hd_data = '../data/6-HD_48350-AD.csv'

    # Read and parse input data
    adx, ady = DataReader(ad_data).get_data()
    hdx, hdy = DataReader(hd_data).get_data()

    def inverse_argmax(l):
        m = 2
        for i in range(len(l)):
            hold = l[i]
            l[i] = [0] * m
            l[i][hold] = 1


    inverse_argmax(ady)
    inverse_argmax(hdy)

    ad_partition = partition_data(adx, ady, [0.8, 0.2])
    ad_batches = CrossValidation(ad_partition[0][0], ad_partition[0][1], 10)
    adx_test, ady_test = np.array(ad_partition[1][0]), np.array(ad_partition[1][1])


    hd_partition = partition_data(hdx, hdy, [0.8, 0.2])
    hd_batches = CrossValidation(hd_partition[0][0], hd_partition[0][1], 10)
    hdx_test, hdy_test = np.array(hd_partition[1][0]), np.array(hd_partition[1][1])

    batch_size = 50 

    n_input = len(adx_test[0])
    # Tensorflow variables

    # Shared layer for MTL hard parameter sharing
    shared_layer_weights = tf.Variable(tf.random_normal([n_input, n_hidden_shared]))
    shared_layer_biases = tf.Variable(tf.random_normal([n_hidden_shared]))

    ad_weights = {
        'hidden_task': tf.Variable(tf.random_normal([n_input, n_hidden_task])),
        'hidden_shared': shared_layer_weights,
        'hidden_2': tf.Variable(tf.random_normal([n_hidden_shared + n_hidden_task, n_hidden_2])),
        'out': tf.Variable(tf.random_normal([n_hidden_2, n_classes]))
    }

    ad_biases = {
        'hidden_task': tf.Variable(tf.random_normal([n_hidden_task])),
        'hidden_shared': shared_layer_biases,
        'hidden_2': tf.Variable(tf.random_normal([n_hidden_2])),
        'out': tf.Variable(tf.random_normal([n_classes]))
    }

    hd_weights = {
        'hidden_task': tf.Variable(tf.random_normal([n_input, n_hidden_task])),
        'hidden_shared': shared_layer_weights,
        'hidden_2': tf.Variable(tf.random_normal([n_hidden_shared + n_hidden_task, n_hidden_2])),
        'out': tf.Variable(tf.random_normal([n_hidden_2, n_classes]))
    }

    hd_biases = {
        'hidden_task': tf.Variable(tf.random_normal([n_hidden_task])),
        'hidden_shared': shared_layer_biases,
        'hidden_2': tf.Variable(tf.random_normal([n_hidden_2])),
        'out': tf.Variable(tf.random_normal([n_classes]))
    }

    # Placeholders
    x = tf.placeholder(tf.float32, [None, n_input])
    y = tf.placeholder(tf.float32, [None, n_classes])
    dim = tf.placeholder(tf.int32)

    model = mtl_2 

    ad_pred = model(x, ad_weights, ad_biases, act=act_func)
    ad_cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=ad_pred, labels=y))


    hd_pred = model(x, hd_weights, hd_biases, act=act_func)
    hd_cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=hd_pred, labels=y))


    global_step = tf.Variable(0, trainable=False)
    starter_learning_rate = 0.1
    learning_rate = tf.train.exponential_decay(starter_learning_rate, 
                                            global_step,
                                            50, 0.5, staircase=False)

    ad_optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(ad_cost)
    hd_optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(hd_cost)

    ad_correct_pred = tf.equal(tf.argmax(ad_pred, 1), tf.argmax(y, 1))
    hd_correct_pred = tf.equal(tf.argmax(hd_pred, 1), tf.argmax(y, 1))
    ad_accuracy = tf.reduce_mean(tf.cast(ad_correct_pred, tf.float32))
    hd_accuracy = tf.reduce_mean(tf.cast(hd_correct_pred, tf.float32))

    init = tf.global_variables_initializer()
    v_acc_hist = []
    max_v_acc = 0
    # Train
    with tf.Session() as sess:
        sess.run(init)
        steps = 2001
        #print('Iteration | Valid loss | Valid acc')
        for i in range(steps):
            if i % 2 == 0:
                train_x, train_y = ad_batches.next_batch(batch_size)
                sess.run(ad_optimizer, feed_dict={x: train_x, y: train_y, global_step: i})
                test_x, test_y = adx_test, ady_test
                v_acc = sess.run(ad_accuracy, feed_dict={x: test_x, y: test_y, dim: len(test_x)})
                v_lss = sess.run(ad_cost, feed_dict={x: test_x, y: test_y, dim: len(test_x)})
                #print('{0:9} |   {1:2.5f} |   {2:2.5f}'.format(i, v_lss, v_acc))
                v_acc_hist.append(v_acc)
            else:
                train_x, train_y = hd_batches.next_batch(batch_size)
                sess.run(hd_optimizer, feed_dict={x: train_x, y: train_y, global_step: i})
                test_x, test_y = hdx_test, hdy_test
                v_acc = sess.run(hd_accuracy, feed_dict={x: test_x, y: test_y, dim: len(test_x)})
                v_lss = sess.run(hd_cost, feed_dict={x: test_x, y: test_y, dim: len(test_x)})
                #print('{0:9} |   {1:2.5f} |   {2:2.5f}'.format(i, v_lss, v_acc))

    print(params, sum(v_acc_hist[-10:]) / 10)
    # Plot validation set and test set loss over epochs
    #x_r = range(0, 1001)
    #fig = plt.figure()
    #ax1 = fig.add_subplot(111)
    #ax1.scatter(x_r, v_lss_hist, label='Validation loss')
    #plt.legend(loc='upper left')
    #plt.show()

shared_search = [128, 256, 512, 1024]
task_vals = [128, 256, 512, 1024]
hidden_2_vals = [128, 256]
act_func_vals = [tf.tanh, tf.nn.relu, tf.sigmoid]
'''
for a in act_func_vals:
    for s in shared_search:
        for t in task_vals:
            for h in hidden_2_vals:
                gridsearch({
                    'n_hidden_shared': s,
                    'n_hidden_task': t,
                    'n_hidden_2': h,
                    'act_func': a
                })
'''
gridsearch({
    'n_hidden_shared': 256,
    'n_hidden_task': 256,
    'n_hidden_2': 256,
    'act_func': tf.tanh
})

