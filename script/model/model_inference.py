# Author: [ChengPeng]
import torch
import pandas as pd
from torch.utils.data import DataLoader, TensorDataset
import torch.nn.functional as F
import numpy as np


def one_hot_encode(sequences, max_len=201):
    encoding = {
        "A": [1, 0, 0, 0],
        "C": [0, 1, 0, 0],
        "G": [0, 0, 1, 0],
        "T": [0, 0, 0, 1],

        "R": [0.5, 0, 0.5, 0],  # A or G
        "Y": [0, 0.5, 0, 0.5],  # C or T
        "S": [0, 0.5, 0.5, 0],  # G or C
        "W": [0.5, 0, 0, 0.5],  # A or T
        "K": [0, 0, 0.5, 0.5],  # G or T
        "M": [0.5, 0.5, 0, 0],  # A or C

        "B": [0, 1 / 3, 1 / 3, 1 / 3],  # C/G/T
        "D": [1 / 3, 0, 1 / 3, 1 / 3],  # A/G/T
        "H": [1 / 3, 1 / 3, 0, 1 / 3],  # A/C/T
        "V": [1 / 3, 1 / 3, 1 / 3, 0],  # A/C/G
        "N": [0.25, 0.25, 0.25, 0.25]
    }
    encoding = {**encoding, **{k.lower(): v for k, v in encoding.items()}}
    encoded_sequences = []
    for sequence in sequences:
        encoded_seq = [encoding[base] for base in sequence]
        if len(sequence) < max_len:
            encoded_seq += [[0.25, 0.25, 0.25, 0.25]] * (max_len - len(sequence))
        encoded_sequences.append(encoded_seq)
    return np.array(encoded_sequences)

def Data_prepare_val(df, max_len=201):
    sequences = list(df["sequence"].values)
    # Find the maximum sequence length
    encoded_sequences = one_hot_encode(sequences, max_len)
    X_val = torch.tensor(encoded_sequences, dtype=torch.float32).unsqueeze(1)
    val_dataset = TensorDataset(X_val)
    val_loader = DataLoader(val_dataset, batch_size=128, shuffle=False)
    return val_loader


def val_model(model, val_data, device="cpu", threshold=0.5):
    model.eval()
    predicted_probs = []
    predicted_labels_list = []
    with torch.inference_mode():
        for batch_idx, (X,) in enumerate(val_data):
            X = X.to(device)
            y_pred = model(X)
            probs = F.softmax(y_pred, dim=1)
            pos_probs = probs[:, 0].cpu().numpy()
            predicted_probs.extend(pos_probs)
            predicted_labels = (pos_probs <= threshold).astype(int)
            predicted_labels_list.extend(predicted_labels)
    return predicted_probs, predicted_labels_list


def filter_sequences_with_model(sequences, model, max_len=201):
    df = pd.DataFrame(sequences, columns=["chrom", "start", "end", "gene", "strand", "sequence"])
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    val_loader = Data_prepare_val(df,max_len)
    prob, labels = val_model(model, val_loader)
    df["label"] = labels  
    df["sequence"] = 1
    filtered_sequences = df[df["label"] == 0][["chrom", "start", "end","sequence" ,"gene", "strand"]]
    filtered_sequences.reset_index(drop=True, inplace=True)
    return filtered_sequences



