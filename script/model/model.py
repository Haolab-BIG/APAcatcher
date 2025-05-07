# Author: [ChengPeng]
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
import torch.nn.functional as F
from einops import rearrange

class PAS_CNN(nn.Module):
    def __init__(self):
        super(PAS_CNN, self).__init__()
        self.conv1 = nn.Conv2d(in_channels=1, out_channels=96, kernel_size=(6, 4))
        self.pool1 = nn.MaxPool2d(kernel_size=(2, 1))
        self.lstm = nn.LSTM(input_size=96, hidden_size=64, num_layers=1, batch_first=True,bidirectional=True)
        self.relu = nn.ReLU()
        self.fc1 = nn.Linear(6272*2, 256)  #
        self.fc2 = nn.Linear(256, 2)

    def forward(self, x):
        x = self.relu(self.conv1(x))
        x = self.pool1(x)
        x = rearrange(x, 'b c h w -> b h (w c)')
        x, _ = self.lstm(x)
        x = x.reshape(x.size(0), -1)
        x = self.relu(self.fc1(x))
        x = self.fc2(x)
        return x

