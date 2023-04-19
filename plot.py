import sys
import pandas as pd
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(
    prog="IsingPlotter",
    description="Ploting of results for Ising Model MC simulation",
)

parser.add_argument('-f', '--file')
parser.add_argument('-o', '--output')

args = parser.parse_args()

inp = args.file or sys.stdin
out = args.output or "graph"

df = pd.read_csv(inp)
fig, axs = plt.subplots(2, 2, figsize=(16, 12))
axs[0, 0].plot(df["temperature"], df["energy"])
axs[0, 0].title.set_text("Energy")
axs[0, 1].plot(df["temperature"], df["specific_heat"])
axs[0, 1].title.set_text("Specific Heat")
axs[1, 0].plot(df["temperature"], df["magnetisation"])
axs[1, 0].title.set_text("Magnetisation")
axs[1, 1].plot(df["temperature"], df["susceptibility"])
axs[1, 1].title.set_text("Susceptibility")
plt.subplots_adjust(wspace=0.1, hspace=0.2)
plt.savefig(out)