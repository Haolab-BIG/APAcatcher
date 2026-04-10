import torch
import pandas as pd
from torch.utils.data import DataLoader, TensorDataset
import torch.nn.functional as F
import numpy as np


def one_hot_encode(sequences, max_len=201):
    lookup = np.full((256, 4), 0.25, dtype=np.float32)
    lookup[ord('A')] = [1, 0, 0, 0]
    lookup[ord('C')] = [0, 1, 0, 0]
    lookup[ord('G')] = [0, 0, 1, 0]
    lookup[ord('T')] = [0, 0, 0, 1]
    lookup[ord('R')] = [0.5, 0, 0.5, 0]
    lookup[ord('Y')] = [0, 0.5, 0, 0.5]
    lookup[ord('S')] = [0, 0.5, 0.5, 0]
    lookup[ord('W')] = [0.5, 0, 0, 0.5]
    lookup[ord('K')] = [0, 0, 0.5, 0.5]
    lookup[ord('M')] = [0.5, 0.5, 0, 0]
    
    for char in "ACGTRYSWKM":
        lookup[ord(char.lower())] = lookup[ord(char)]

    num_samples = len(sequences)
    encoded_array = np.zeros((num_samples, max_len, 4), dtype=np.float32)
    
    for i, seq in enumerate(sequences):
        curr_seq = seq[:max_len]
        curr_len = len(curr_seq)
        encoded_array[i, :curr_len, :] = lookup[np.frombuffer(curr_seq.encode('ascii'), dtype=np.uint8)]
        encoded_array[i, curr_len:, :] = 0.25
        
    return encoded_array


def Data_prepare_val(sequences, max_len=201):
    encoded_sequences = one_hot_encode(sequences, max_len)
    X_val = torch.from_numpy(encoded_sequences).float().unsqueeze(1)
    val_dataset = TensorDataset(X_val)
    return DataLoader(val_dataset, batch_size=256, shuffle=False, num_workers=4)



def val_model(model, val_loader, device="cpu", threshold=0.5):
    model.to(device)
    model.eval()
    all_pos_probs = []
    with torch.inference_mode():
        for (X,) in val_loader:
            X = X.to(device)
            logits = model(X)
            probs = F.softmax(logits, dim=1)
            pos_probs = probs[:, 0].cpu().numpy() 
            all_pos_probs.append(pos_probs)
    predicted_probs = np.concatenate(all_pos_probs)
    predicted_labels = (predicted_probs <= threshold).astype(int)
    return predicted_probs, predicted_labels


def filter_sequences_with_model(sequences_info, model, max_len=201, device="cpu"):
    raw_sequences = [item[5] for item in sequences_info] 
    val_loader = Data_prepare_val(raw_sequences, max_len)
    _, labels = val_model(model, val_loader, device=device)
    df = pd.DataFrame(sequences_info, columns=["chrom", "start", "end", "gene", "strand","sequence"])
    df["label"] = labels
    df["score"] = 1
    filtered_df = df[df["label"] == 0].copy()
    result = filtered_df[["chrom", "start", "end", "score", "gene", "strand"]]
    return result


