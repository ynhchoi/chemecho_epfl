from jcamp import jcamp_read
import pandas as pd
from io import StringIO
import nistchempy as nist
import matplotlib.pyplot as plt


def from_df_to_csv(df_spectrum) :
    buffer_csv = StringIO()
    return df_spectrum.to_csv(buffer_csv, index=False)

