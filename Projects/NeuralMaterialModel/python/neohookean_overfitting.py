
import os
from functools import cmp_to_key
from statistics import mode
# os.environ["CUDA_VISIBLE_DEVICES"] = "0"
os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = "true"
import math
import numpy as np
import tensorflow as tf
from model import *
from tensorflow.keras.optimizers import Adam
from tensorflow.keras import Model
from tensorflow.keras.models import load_model
import matplotlib.pyplot as plt

tf.keras.backend.set_floatx('float64')

def relativeL2(y_true, y_pred):
    # msle = tf.keras.losses.MeanSquaredLogarithmicError()
    if (y_true.shape[1] > 1):
        stress_norm = tf.norm(y_true, ord='euclidean', axis=1)
        norm = tf.tile(tf.keras.backend.expand_dims(stress_norm, 1), tf.constant([1, 3]))
        y_true_normalized = tf.divide(y_true, norm)
        y_pred_normalized = tf.divide(y_pred, norm)
        return K.mean(K.square(y_true_normalized - y_pred_normalized))
        # exit(0)
    else:
        
        y_true_normalized = tf.ones(y_true.shape, dtype=tf.float64)
        y_pred_normalized = tf.divide(y_pred + K.epsilon(), y_true + K.epsilon())
        return K.mean(K.square((y_true_normalized - y_pred_normalized) * tf.constant(1.0, dtype=tf.float64)))

def loadDataSplitTest(filename, shuffle = True, ignore_unconverging_result = True):
    all_data = []
    all_label = [] 
    
    for line in open(filename).readlines():
        item = [float(i) for i in line.strip().split(" ")]
        if (ignore_unconverging_result):
            if (item[-1] > 1e-6 or math.isnan(item[-1])):
                continue
        data = item[:2]
        data.append(2.0 * item[2])

        label = item[3:7]
        
        all_data.append(data)
        all_label.append(label)
        
    
    all_data = np.array(all_data[:]).astype(np.float64)
    all_label = np.array(all_label[:]).astype(np.float64)
    indices = np.arange(all_data.shape[0])
    if (shuffle):
        np.random.shuffle(indices)
    all_data = all_data[indices]
    all_label = all_label[indices]
    

    return all_data, all_label


def generator(train_data, train_label):    
    indices = np.arange(train_data.shape[0])
    while True:
        np.random.shuffle(indices)
        yield train_data[indices], train_label[indices]

w_grad = tf.constant(1.0, dtype=tf.float64)
w_e = tf.constant(1.0, dtype=tf.float64)

@tf.function
def trainStep(opt, lambdas, sigmas, model, train_vars):
    
    
    with tf.GradientTape(persistent=True) as tape:
        tape.watch(train_vars)
        tape.watch(lambdas)
        
        psi = model(lambdas)
        dedlambda = tape.gradient(psi, lambdas)
        batch_dim = psi.shape[0]
        stress_gt = tf.slice(sigmas, [0, 0], [batch_dim, 3])
        potential_gt = tf.slice(sigmas, [0, sigmas.shape[1]-1], [batch_dim, 1])
        stress_pred = tf.slice(dedlambda, [0, 0], [batch_dim, 3])
        
        grad_loss = w_grad * relativeL2(stress_gt, stress_pred)
        e_loss = w_e * relativeL2(potential_gt, psi)

        loss = grad_loss + e_loss
        
    dLdw = tape.gradient(loss, train_vars)
    opt.apply_gradients(zip(dLdw, train_vars))
    gradNorm = tf.math.sqrt(tf.reduce_sum([tf.reduce_sum(gi*gi) for gi in dLdw]))
    # gradNorm = -1.0
    
    del tape
    return grad_loss, e_loss, gradNorm

@tf.function
def testStep(lambdas, sigmas, model):
    
    with tf.GradientTape() as tape:
        tape.watch(lambdas)
        psi = model(lambdas)
        dedlambda = tape.gradient(psi, lambdas)
        batch_dim = psi.shape[0]
        stress_gt = tf.slice(sigmas, [0, 0], [batch_dim, 3])
        potential_gt = tf.slice(sigmas, [0, sigmas.shape[1]-1], [batch_dim, 1])
        stress_pred = tf.slice(dedlambda, [0, 0], [batch_dim, 3])
        
        grad_loss = w_grad * relativeL2(stress_gt, stress_pred)
        e_loss = w_e * relativeL2(potential_gt, psi)
    del tape
    return grad_loss, e_loss, stress_pred, psi

@tf.function
def valueGradHessian(inputs, model):
    batch_dim = inputs.shape[0]
    with tf.GradientTape() as tape_outer:
        tape_outer.watch(inputs)
        with tf.GradientTape() as tape:
            tape.watch(inputs)
            psi = model(inputs, training=False)
            dedlambda = tape.gradient(psi, inputs)
            
            stress = tf.slice(dedlambda, [0, 0], [batch_dim, 3])
            
    C = tape_outer.batch_jacobian(stress, inputs)[:, :, :]
    del tape
    del tape_outer
    return psi, stress, C

def NHAutodiffTest(inputs, lambda_tf, mu_tf, data_type = tf.float64):
    strain = inputs
    batch_dim = strain.shape[0]
    with tf.GradientTape() as tape_outer:
        tape_outer.watch(strain)
        with tf.GradientTape() as tape:
            tape.watch(strain)
            strain_xx = tf.gather(strain, [0], axis = 1)
            strain_yy = tf.gather(strain, [1], axis = 1)
            
            strain_xy = tf.constant(0.5, dtype=data_type) * tf.gather(strain, [2], axis = 1)
            strain_vec_reorder = tf.concat((strain_xx, strain_xy, strain_xy, strain_yy), axis=1)
            
            strain_tensor = tf.reshape(strain_vec_reorder, (batch_dim, 2, 2))    
                        
            righCauchy = tf.constant(2.0, dtype=data_type) * strain_tensor + tf.eye(2, batch_shape=[batch_dim], dtype=data_type)
            
            J = tf.math.sqrt(tf.linalg.det(righCauchy))
            
            I1 = tf.linalg.trace(righCauchy)
            C1 = tf.constant(0.5 * mu_tf, dtype=data_type)
            D1 = tf.constant(lambda_tf * 0.5, dtype=data_type)
            lnJ = tf.math.log(J)
            psi = C1 * (I1 - tf.constant(2.0, dtype=data_type) - tf.constant(2.0, dtype=data_type) * lnJ) + D1 * (lnJ*lnJ)
            
            stress = tape.gradient(psi, strain)
            # print(stress)
            # exit(0)
    C = tape_outer.batch_jacobian(stress, strain)
    print(C)
    exit(0)

@tf.function
def psiGradHessNH(strain, data_type = tf.float64):
    lambda_tf = 26.0 * 0.48 / (1.0 + 0.48) / (1.0 - 2.0 * 0.48)
    mu_tf = 26.0 / 2.0 / (1.0 + 0.48)
    youngsModulus = 26.0
    poissonsRatio = 0.48

    batch_dim = strain.shape[0]
    with tf.GradientTape() as tape_outer:
        tape_outer.watch(strain)
        with tf.GradientTape() as tape:
            tape.watch(strain)
            
            strain_xx = tf.gather(strain, [0], axis = 1)
            strain_yy = tf.gather(strain, [1], axis = 1)
            
            strain_xy = tf.constant(0.5, dtype=data_type) * tf.gather(strain, [2], axis = 1)
            strain_vec_reorder = tf.concat((strain_xx, strain_xy, strain_xy, strain_yy), axis=1)
            
            strain_tensor = tf.reshape(strain_vec_reorder, (batch_dim, 2, 2))    
                        
            righCauchy = tf.constant(2.0, dtype=data_type) * strain_tensor + tf.eye(2, batch_shape=[batch_dim], dtype=data_type)
            
            J = tf.math.sqrt(tf.linalg.det(righCauchy))
            
            I1 = tf.linalg.trace(righCauchy)
            C1 = tf.constant(0.5 * mu_tf, dtype=data_type)
            D1 = tf.constant(lambda_tf * 0.5, dtype=data_type)
            lnJ = tf.math.log(J)
            psi = C1 * (I1 - tf.constant(2.0, dtype=data_type) - tf.constant(2.0, dtype=data_type) * lnJ) + D1 * (lnJ*lnJ)
            
            stress = tape.gradient(psi, strain)
            # print(stress)
            # exit(0)
    C = tape_outer.batch_jacobian(stress, strain)
    return psi, stress, C

def psiGradHessStVK(strain, data_type = tf.float64):
    lambda_tf = 26.0 * 0.48 / (1.0 + 0.48) / (1.0 - 2.0 * 0.48)
    mu_tf = 26.0 / 2.0 / (1.0 + 0.48)
    youngsModulus = 26.0
    poissonsRatio = 0.48

    # NHAutodiffTest(tf.gather(strain, [5], axis=0), lambda_tf, mu_tf)

    f = youngsModulus / (1 - poissonsRatio * poissonsRatio)
    stiffnessTensor = np.reshape(np.array([f, f * poissonsRatio, 0, f * poissonsRatio, f, 0, 0, 0, (1 - poissonsRatio) / 2 * f]), (3,3))
    print("isotropic C ", stiffnessTensor)
	

    batch_dim = strain.shape[0]
    with tf.GradientTape() as tape_outer:
        tape_outer.watch(strain)
        with tf.GradientTape() as tape:
            tape.watch(strain)
            
            strain_xx = tf.gather(strain, [0], axis = 1)
            strain_yy = tf.gather(strain, [1], axis = 1)
            
            strain_xy = tf.constant(0.5, dtype=data_type) * tf.gather(strain, [2], axis = 1)
            strain_vec_reorder = tf.concat((strain_xx, strain_xy, strain_xy, strain_yy), axis=1)
            
            strain_tensor = tf.reshape(strain_vec_reorder, (batch_dim, 2, 2))    
            
            E2 = tf.matmul(strain_tensor, strain_tensor)
            psi = tf.constant(0.5, dtype=data_type) *tf.math.pow(tf.linalg.trace(strain_tensor), tf.constant(2.0, dtype=data_type)) + tf.linalg.trace(E2)

            stress = tape.gradient(psi, strain)
            # print(stress)
            # exit(0)
    C = tape_outer.batch_jacobian(stress, strain)
    return C, stress, psi

@tf.function
def computeDirectionalStiffnessNH(inputs, thetas):
    batch_dim = inputs.shape[0]
    thetas = tf.expand_dims(thetas, axis=1)

    d_voigt = tf.concat((tf.math.cos(thetas) * tf.math.cos(thetas), 
                        tf.math.sin(thetas) * tf.math.sin(thetas), 
                        tf.math.sin(thetas) * tf.math.cos(thetas)), 
                        axis = 1)

    
    psi, stress, C = psiGradHessNH(tf.convert_to_tensor(inputs))
   
    Sd = tf.linalg.matvec(tf.linalg.inv(C[0, :, :]), d_voigt[0, :])
    dTSd = tf.expand_dims(tf.tensordot(d_voigt[0, :], Sd, 1), axis=0)
    
    for i in range(1, C.shape[0]):
        
        Sd = tf.linalg.matvec(tf.linalg.inv(C[i, :, :]), d_voigt[i, :])
        dTSd = tf.concat((tf.expand_dims(tf.tensordot(d_voigt[i, :], Sd, 1), axis=0), dTSd), 0)
        
    stiffness = tf.squeeze(tf.math.divide(tf.ones((batch_dim), dtype=tf.float64), tf.expand_dims(dTSd, axis=0)))
    
    return stiffness


@tf.function
def computeDirectionalStiffness(inputs, thetas, model):
    batch_dim = inputs.shape[0]
    thetas = tf.expand_dims(thetas, axis=1)

    d_voigt = tf.concat((tf.math.cos(thetas) * tf.math.cos(thetas), 
                        tf.math.sin(thetas) * tf.math.sin(thetas), 
                        tf.math.sin(thetas) * tf.math.cos(thetas)), 
                        axis = 1)

    
    psi, stress, C = valueGradHessian(inputs, model)
   
    Sd = tf.linalg.matvec(tf.linalg.inv(C[0, :, :]), d_voigt[0, :])
    dTSd = tf.expand_dims(tf.tensordot(d_voigt[0, :], Sd, 1), axis=0)
    
    for i in range(1, C.shape[0]):
        
        Sd = tf.linalg.matvec(tf.linalg.inv(C[i, :, :]), d_voigt[i, :])
        dTSd = tf.concat((tf.expand_dims(tf.tensordot(d_voigt[i, :], Sd, 1), axis=0), dTSd), 0)
        
    stiffness = tf.squeeze(tf.math.divide(tf.ones((batch_dim), dtype=tf.float64), tf.expand_dims(dTSd, axis=0)))
    eng_strain = 0.1
    stiffness2 = tf.constant(2.0, dtype=tf.float64) * tf.math.divide(tf.squeeze(psi), tf.math.pow(tf.constant(eng_strain, dtype=tf.float64), tf.constant(2.0, dtype=tf.float64)) * tf.ones((batch_dim), dtype=tf.float64))
    
    return stiffness, stiffness2

def toPolarData(half):
    full = half
    n_sp_theta = len(half)
    for i in range(n_sp_theta):
        full = np.append(full, full[i])
    full = np.append(full, full[0])
    return full

def train(model_name, train_data, train_label, validation_data, validation_label):
    batch_size = np.minimum(60000, len(train_data))
    # print("batch size: {}".format(batch_size))
    # model = buildSingleFamilyModel(n_tiling_params)
    model = buildConstitutiveModel()    
    train_vars = model.trainable_variables
    opt = Adam(learning_rate=1e-4)
    max_iter = 80000

    val_lambdasTF = tf.convert_to_tensor(validation_data)
    val_sigmasTF = tf.convert_to_tensor(validation_label)

    losses = [[], []]
    current_dir = os.path.dirname(os.path.realpath(__file__))
    # model.load_weights("/home/yueli/Documents/ETH/WuKong/Projects/Tilisng2D/python/Models/67/" + model_name + '.tf')
    count = 0
    with open('counter.txt', 'r') as f:
        count = int(f.read().splitlines()[-1])
    f = open("counter.txt", "w+")
    f.write(str(count+1))
    f.close()
    
    
    save_path = os.path.join(current_dir, 'Models/' + str(count) + "/")
    if not os.path.exists(save_path):
        os.mkdir(save_path)
    g_norm0 = 0
    iter = 0
    
    for iteration in range(max_iter):
        lambdas, sigmas = next(generator(train_data, train_label))
        if batch_size == -1:
            batch = 1
        else:
            batch = int(np.floor(len(lambdas) / batch_size))
        
        train_loss_grad = 0.0
        train_loss_e = 0.0
        g_norm_sum = 0.0
        for i in range(batch):
            mini_bacth_lambdas = lambdas[i * batch_size:(i+1) * batch_size]
            mini_bacth_sigmas = sigmas[i * batch_size:(i+1) * batch_size]

            lambdasTF = tf.convert_to_tensor(mini_bacth_lambdas)
            sigmasTF = tf.convert_to_tensor(mini_bacth_sigmas)
            
            grad, e, g_norm = trainStep(opt, lambdasTF, sigmasTF, model, train_vars)
            
            train_loss_grad += grad
            train_loss_e += e
            g_norm_sum += g_norm
        if (iteration == 0):
            g_norm0 = g_norm_sum
        validation_loss_grad, validation_loss_e, _, _ = testStep(val_lambdasTF, val_sigmasTF, model)
        
        losses[0].append(train_loss_grad + train_loss_e)
        losses[1].append(validation_loss_grad + validation_loss_e)
        print("epoch: {}/{} train_loss_grad: {} train_loss e: {}, validation_loss_grad:{} loss_e:{} |g|: {}, |g_init|: {} ".format(iteration, max_iter, train_loss_grad, train_loss_e, \
                         validation_loss_grad, validation_loss_e, \
                        g_norm_sum, g_norm0))
        if iteration % 10000 == 0:
            model.save_weights(save_path + model_name + '.tf')
    model.save_weights(save_path + model_name + '.tf')
    return count

def validate(count, model_name, validation_data, validation_label):
    current_dir = os.path.dirname(os.path.realpath(__file__))
    save_path = os.path.join(current_dir, 'Models/' + str(count) + "/")
    # model = loadSingleFamilyModel(n_tiling_params)
    model = buildConstitutiveModel()
    model.load_weights(save_path + model_name + '.tf')
    # model.save(save_path + model_name + '.h5')
    grad_loss, e_loss, sigma, energy = testStep(validation_data, validation_label, model)

    print("validation loss grad: {} energy: {}".format(grad_loss, e_loss)) 



if __name__ == "__main__":
    
    data_file = "/media/yueli/Elements/SandwichStructure/Homo/homo_uni_bi_shuffled.txt"
    # data_file = "/home/yueli/Documents/ETH/SandwichStructure/SampleStrain/homo_sample_theta_1.1.txt"
    
    data_all, label_all = loadDataSplitTest(data_file, shuffle=False, ignore_unconverging_result=True)

    five_percent = int(len(data_all) * 0.05)

    train_data =  data_all[:-five_percent]
    train_label =  label_all[:-five_percent]

    validation_data = data_all[-five_percent:]
    validation_label = label_all[-five_percent:]

    model_name = "homo_NH"
    # exp_id = train(model_name, 
    #     train_data, train_label, validation_data, validation_label)

    model = buildConstitutiveModel()
    # model.load_weights("./Models/9/" + model_name + '.tf')
    # model.save("./Models/9/", save_format='tf')


    class NeuralModel(keras.Model):
        def __init__(self):
            super().__init__()
            self.model = model
            self.model.load_weights("./Models/9/" + model_name + '.tf')
        def call(self, strain):
            with tf.GradientTape() as tape2:
                tape2.watch(strain)
                with tf.GradientTape() as tape:
                    tape.watch(strain)
                    energy_density = self.model(strain)
                gradients = tape.gradient(energy_density, strain)
            hessians = tape2.batch_jacobian(gradients, strain)
            del tape
            del tape2
            return {"energy": energy_density, "gradient": gradients, "hessian": hessians}
    neural_model = NeuralModel()
    output = neural_model(tf.constant([[1.0, 1.0, 1.0]]))
    neural_model.save("./Models/9/", save_format='tf')
    
    