import matplotlib.pyplot as plt
from jcamp import jcamp_read
import pandas as pd
from io import StringIO


def ir_graph(dict_spectrum, cas_number) :
    df_spectrum = pd.DataFrame({"Wavenumber":dict_spectrum.x.tolist(), "Transmittance":dict_spectrum.y.tolist()})
    compound = nist.get_compound(cas_number)

    fig, plt = plt.subplots()
    plt.plot(df_spectrum["Wavenumber"], df_spectrum["Transmittance"], color="black")
    plt.set_xlabel("Wavenumber []")
    plt.set_ylabel("Transmittance [%]")
    plt.set_title(f"IR spectrum of {compound.name}")

    plt.grid(visible=False)
    plt.set_xlim(4000, 400)
    plt.set_ylim(0, 100)

    return fig, plt