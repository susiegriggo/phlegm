# imports 
import numpy as np
import pandas as pd 
from sklearn.model_selection import KFold
from loguru import logger
import glob
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
import joblib
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.optimizers import Adam
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset

logger.add("trainer.log")

def get_af2_emb(af2_embeddings_dir: str, structure_id: int, use_pairwise: bool, id_: str = None) -> np.array:
    """
    Get AF2 embeddings from ColabFold output directory.
    *** CODE FROM https://github.com/labstructbioinf/dc2_oligo *** 

    Parameters:
        af2_embeddings_dir (str): Path to the ColabFold output directory.
        model_id (int): Model ID to retrieve embeddings from.
        use_pairwise (bool): Whether to include pairwise embeddings.
        id_ (str, optional): ID of the protein from dataframe to retrieve embeddings

    Returns:
        np.ndarray: Array containing the AF2 embeddings.
    """
    if id_ is None:
        single_repr_fns = sorted(glob.glob(f"{af2_embeddings_dir}/*_single_repr_rank_*_model_{structure_id+1}_*"))
        pair_repr_fns = sorted(glob.glob(f"{af2_embeddings_dir}/*_pair_repr_rank_*_model_{structure_id+1}_*"))
    else:
        single_repr_fns = sorted(glob.glob(f"{af2_embeddings_dir}/{id_}/*_single_repr_rank_*_model_{structure_id+1}_*"))
        pair_repr_fns = sorted(glob.glob(f"{af2_embeddings_dir}/{id_}/*_pair_repr_rank_*_model_{structure_id+1}_*"))

    mat = np.load(single_repr_fns[0]).mean(axis=0)

    if use_pairwise:
        mat = np.hstack((mat, np.load(pair_repr_fns[0]).mean(axis=0).mean(axis=0)))

    return mat


def logistic_regression(X_train, y_train, max_iter=1000, balanced=True):
    """
    Train a logistic regression model.

    Parameters:
        X_train (np.array): The training data.
        y_train (np.array): The training labels.
        max_iter (int): Maximum number of iterations for the solver.
        balanced (bool): Whether to use class balancing.

    Returns:
        LogisticRegression: The trained logistic regression model.
    """
    class_weight = 'balanced' if balanced else None
    model = LogisticRegression(max_iter=max_iter, random_state=42, class_weight=class_weight)
    model.fit(X_train, y_train)
    return model

def random_forest(X_train, y_train, n_estimators=100, max_depth=None, random_state=42):
    """
    Train a random forest classifier.

    Parameters:
        X_train (np.array): The training data.
        y_train (np.array): The training labels.
        n_estimators (int): The number of trees in the forest.
        max_depth (int, optional): The maximum depth of the tree.
        random_state (int): Random seed for reproducibility.

    Returns:
        RandomForestClassifier: The trained random forest model.
    """
    model = RandomForestClassifier(n_estimators=n_estimators, max_depth=max_depth, random_state=random_state)
    model.fit(X_train, y_train)
    return model

def support_vector_machine(X_train, y_train, C=1.0, kernel='linear', random_state=42):
    """
    Train a support vector machine classifier.

    Parameters:
        X_train (np.array): The training data.
        y_train (np.array): The training labels.
        C (float): Regularization parameter.
        kernel (str): Specifies the kernel type to be used in the algorithm.
        random_state (int): Random seed for reproducibility.

    Returns:
        SVC: The trained support vector machine model.
    """
    model = SVC(C=C, kernel=kernel, probability=True, random_state=random_state)
    model.fit(X_train, y_train)
    return model

def shallow_neural_network(X_train, y_train, input_dim, num_classes, epochs=50, batch_size=32, learning_rate=0.001):
    """
    Train a shallow neural network using Keras.

    Parameters:
        X_train (np.array): The training data.
        y_train (np.array): The training labels.
        input_dim (int): The number of input features.
        num_classes (int): The number of output classes.
        epochs (int): The number of epochs to train the model.
        batch_size (int): The number of samples per gradient update.
        learning_rate (float): Learning rate for the optimizer.

    Returns:
        Sequential: The trained Keras model.
    """
    model = Sequential()
    model.add(Dense(64, input_dim=input_dim, activation='relu'))
    model.add(Dense(32, activation='relu'))
    model.add(Dense(num_classes, activation='softmax'))
    
    model.compile(optimizer=Adam(learning_rate=learning_rate), loss='sparse_categorical_crossentropy', metrics=['accuracy'])
    model.fit(X_train, y_train, epochs=epochs, batch_size=batch_size, verbose=1)
    
    return model

class ShallowNN(nn.Module):
    """
    A shallow neural network using PyTorch.

    Parameters:
        input_dim (int): The number of input features.
        num_classes (int): The number of output classes.
    """
    def __init__(self, input_dim, num_classes):
        super(ShallowNN, self).__init__()
        self.fc1 = nn.Linear(input_dim, 64)
        self.fc2 = nn.Linear(64, 32)
        self.fc3 = nn.Linear(32, num_classes)

    def forward(self, x):
        x = torch.relu(self.fc1(x))
        x = torch.relu(self.fc2(x))
        x = torch.softmax(self.fc3(x), dim=1)
        return x

def shallow_neural_network_pytorch(X_train, y_train, input_dim, num_classes, epochs=50, batch_size=32, learning_rate=0.001):
    """
    Train a shallow neural network using PyTorch.

    Parameters:
        X_train (np.array): The training data.
        y_train (np.array): The training labels.
        input_dim (int): The number of input features.
        num_classes (int): The number of output classes.
        epochs (int): The number of epochs to train the model.
        batch_size (int): The number of samples per gradient update.
        learning_rate (float): Learning rate for the optimizer.

    Returns:
        ShallowNN: The trained PyTorch model.
    """
    model = ShallowNN(input_dim, num_classes)
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)

    train_dataset = TensorDataset(torch.tensor(X_train, dtype=torch.float32), torch.tensor(y_train, dtype=torch.long))
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)

    for epoch in range(epochs):
        model.train()
        for X_batch, y_batch in train_loader:
            optimizer.zero_grad()
            outputs = model(X_batch)
            loss = criterion(outputs, y_batch)
            loss.backward()
            optimizer.step()

    return model

def train_model(af2_embeddings, symmetries, k_folds, num_structures, use_pairwise=True, output_dir=None, model_type='logistic_regression'): 
    # train the model 
    # af2_embeddings: the embeddings of the af2 
    # symmetries: the symmetries of the af2 
    # k_folds: the number of folds for cross validation
    # num_models: the number of structures included in the alphafold embeddings
    # return the model 

    # load in the symmmetries: 
    sym_df = pd.read_csv(symmetries, '\t')
    num_classes = len(sym_df['symmetry'].unique())

    # initalse arrays to store the results
    results = np.zeros((num_structures, len(sym_df), num_classes))
    internal_representations = np.zeros((num_structures, len(sym_df)))
    model = {}
    fold_indices = {}

    for s in range(num_structures):  # loop through the number of structures that are present - there are 5 structures for each embedding currently

        logger.info(f"Training model for structure {s} of {num_structures}")

        # get the embeddings 
        X = np.asarray([get_af2_emb(af2_embeddings, structure_id=s, use_pairwise=use_pairwise, id_=id_) for id_ in sym_df['pdb'].values])
        y = sym_df['symmetry'].values

        # intialise Kfold cross-validation 
        cv = KFold(n_splits=k_folds, shuffle=True, random_state=42)

        # loop through this fold 
        for k, (train_idx, test_idx) in enumerate(cv.split(X, y)): 

            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]

            if model_type == 'logistic_regression':
                model = logistic_regression(X_train, y_train)
            elif model_type == 'random_forest':
                model = random_forest(X_train, y_train)
            elif model_type == 'svm':
                model = support_vector_machine(X_train, y_train)
            elif model_type == 'shallow_nn':
                model = shallow_neural_network_pytorch(X_train, y_train, input_dim=X_train.shape[1], num_classes=num_classes)
            else:
                raise ValueError(f"Unsupported model type: {model_type}")

            # store the predicted probabilities and internal representations
            if model_type == 'shallow_nn':
                model.eval()
                with torch.no_grad():
                    results[s, test_idx, :] = model(torch.tensor(X_test, dtype=torch.float32)).numpy()
                    internal_representations[s, test_idx] = model(torch.tensor(X_test, dtype=torch.float32)).argmax(axis=1).numpy()
            else:
                results[s, test_idx, :] = model.predict_proba(X_test)
                internal_representations[s, test_idx] = model.decision_function(X_test)
            model[f"clf_{s}_{k}"] = model

            # Save fold indices
            fold_indices[f"structure_{s}_fold_{k}"] = {'train_idx': train_idx, 'test_idx': test_idx}

            # Compute metrics to evaluate the quality of the model's performance
            y_pred = model.predict(X_test) if model_type != 'shallow_nn' else model(torch.tensor(X_test, dtype=torch.float32)).argmax(axis=1).numpy()
            accuracy = accuracy_score(y_test, y_pred)
            precision = precision_score(y_test, y_pred, average='weighted')
            recall = recall_score(y_test, y_pred, average='weighted')
            f1 = f1_score(y_test, y_pred, average='weighted')

            logger.info(f'Structure {s}, Fold {k} - Accuracy: {accuracy}') 
            logger.info(f'Structure {s}, Fold {k} - Precision: {precision}')
            logger.info(f'Structure {s}, Fold {k} - Recall: {recall}')
            logger.info(f'Structure {s}, Fold {k} - F1 Score: {f1}')

    # save to a file 
    if output_dir: 
        joblib.dump(model, f"{output_dir}/model_{s}_{k}.joblib")
        np.save(f'{output_dir}/internal_representations_{s}_{k}.npy', internal_representations[s, test_idx]) #need to be cautious of only using predictions out of fold 
        np.save(f'{output_dir}/fold_indices.npy', fold_indices)

    return model, results, internal_representations









