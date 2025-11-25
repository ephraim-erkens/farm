##################################################################################################
### collection of functions that build and fit neural networks

# example use:
# import neural_networks as nn
# ffc = nn.FeedForwardClassification
# model, hist = ffc(...)

# when changes are applied to neural_networks.py, the module needs to be reloaded, e.g. using:
# import importlib
# importlib.reload(neural_networks)


##################################################################################################

import pandas as pd
import numpy as np
import time

from sklearn.metrics import confusion_matrix

import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Activation, Dropout
from tensorflow.keras import regularizers
from tensorflow.keras.optimizers import Adam, SGD

from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import train_test_split


class FNNClassifier:

    
    def __init__(self, input_shape, compound_name, num_classes=5, class_names=[1,2,3,4,5]):
        self.input_shape = input_shape
        self.compound_name = compound_name
        self.num_classes = num_classes
        self.class_names = class_names
        self.name = 'Feedforward Neural Network'
        self.short_name='FNN'
        self.model = None
        self.history = None
        self.confusion = None
        self.fit_time = None
        self.performance = None

    
    def build_model(
        self,
        optimizer='adam',
        # learning_rate=None,
        activation='sigmoid',
        layers=[128, 128, 64],
        regularizers=[None, None, regularizers.l2(0.01)] # we may want to use l1 or elastic net as default
    ):


        # N_VALIDATION = X_val.shape[0]
        # N_TRAIN = X_train.shape[0]
        # BUFFER_SIZE = int(1e4)
        # BATCH_SIZE = 500 # modify?
        STEPS_PER_EPOCH = 7 # N_TRAIN//BATCH_SIZE
        lr_schedule = tf.keras.optimizers.schedules.InverseTimeDecay(
            0.001,
            decay_steps=STEPS_PER_EPOCH*1000,
            decay_rate=1,
            staircase=False
        )
        
        self.model = Sequential()
        self.model.add(tf.keras.Input(shape=self.input_shape))

        if len(layers) != len(regularizers):
            raise ValueError('layers and regularizers lists must have same lenghts')
        for num_nodes, regularizer in zip(layers, regularizers):
            self.model.add(Dense(num_nodes, activation=activation, kernel_regularizer=regularizer))
            self.model.add(Dropout(0.5))
        if self.num_classes==2:
            self.model.add(Dense(1, activation='sigmoid')) ## final prediction layer
        elif self.num_classes>2:
            self.model.add(Dense(self.num_classes, activation='softmax')) ## final prediction layer
        
        if optimizer.lower()=='adam':
            optimizer = Adam(learning_rate=lr_schedule)
        else:
            raise ValueError('Unknown optimizer')

        if self.num_classes==2:
            self.model.compile(optimizer=optimizer, loss='binary_crossentropy', metrics=['accuracy', 'precision', 'recall', 'AUC'])
        elif self.num_classes>2:
            self.model.compile(optimizer=optimizer, loss='categorical_crossentropy', metrics=['accuracy', 'precision', 'recall', 'f1_score'])


    def train(
        self,
        X_train,
        y_train,
        X_val,
        y_val,
        class_weight=None,
        max_epochs=500,
        batch_size=32,
        patience=50,
        verbose=0,
    ):

        callbacks = [
            tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=patience, restore_best_weights=True)
        ]

        start = time.time()
        hist = self.model.fit(
            X_train, y_train,
            epochs=max_epochs, 
            batch_size=batch_size, 
            validation_data=(X_val, y_val), 
            callbacks=callbacks, 
            verbose=verbose,
            class_weight=class_weight,
        )
        end = time.time()
        self.fit_time = end-start
        self.history = pd.DataFrame(hist.history)


    def evaluate(
        self,
        X_test,
        y_test,
        verbose=0,
        threshold=0.5
    ):

        y_predicted = self.model.predict(X_test, verbose=verbose)
        if self.num_classes==2:
            y_predicted = np.array([1 if x >= threshold else 0 for x in y_predicted]) # threshhold may be set differently here after ROC review
            confusion = confusion_matrix(y_test, y_predicted)
        elif self.num_classes>2:
            y_test_indices = np.argmax(y_test, axis=1)
            y_pred_indices = np.argmax(y_predicted, axis=1)
            confusion = confusion_matrix(y_test_indices, y_pred_indices)
        self.confusion = confusion

        metrics = self.model.evaluate(X_test, y_test, return_dict=True, verbose=0)
        
        # store performances
        log_dict = dict()
        log_dict['compound'] = self.compound_name
        log_dict['train_accuracy'] = self.history.accuracy.iloc[-1]
        log_dict['val_accuracy'] = self.history.val_accuracy.iloc[-1]
        log_dict['test_accuracy'] = metrics['accuracy']
        log_dict['train_precision'] = self.history.precision.iloc[-1]
        log_dict['val_precision'] = self.history.val_precision.iloc[-1]
        log_dict['test_precision'] = metrics['precision']
        log_dict['train_recall'] = self.history.recall.iloc[-1]
        log_dict['val_recall'] = self.history.val_recall.iloc[-1]
        log_dict['test_recall'] = metrics['recall']
        if self.num_classes==2:
            log_dict['train_auc'] = self.history.AUC.iloc[-1]
            log_dict['val_auc'] = self.history.val_AUC.iloc[-1]
            log_dict['test_auc'] = metrics['AUC']
        elif self.num_classes>2:
            for j in range(self.num_classes):
                log_dict[f'test_f1_class{j+1}'] = float(self.history.f1_score.iloc[-1][j])
        for row in range(confusion.shape[0]):
            for col in range(confusion.shape[1]):
                log_dict[f'cm{row}{col}'] = confusion[row][col]

        self.performance = log_dict
    
        return log_dict
        

