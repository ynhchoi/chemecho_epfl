import streamlit as st
import nistchempy as nist

compound = nist.get_compound("67-64-1")
print(compound)