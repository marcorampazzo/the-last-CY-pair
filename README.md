# The last CY pair: computations for the paper

This repository contains the Python and Macaulay2 computations accompanying *The Last Picard Rank 1 Double-Mirror Calabi-Yau Pair?* by Michał Kapustka, Marco Rampazzo, and Prajwal Samal, and describes how to run the Python code. 

### What's inside

#### Macaulay2
- `mirror_principal_period.m2` refers to Appendix A2
- `mirror_PF_operator.m2` refers to Appendix A3

To install Macaulay2 (M2), follow the [official installation instructions](https://github.com/Macaulay2/M2/wiki) for your operating system: macOS, Linux, or Windows (via Windows Subsystem for Linux). After installation, type `M2` in a terminal to start it. These scripts do not require the authors' Python library; for help using Macaulay2, see its [documentation](https://macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/Macaulay2Doc/html/).

#### Python
- `hodge_numbers.py` refers to Appendix A1
- `pretilting.py` refers to Appendix A4
- `open_cohomology.py` refers to Appendix A5

These scripts use the authors' [homogeneous-varieties library](https://github.com/marcorampazzo/homogeneous-varieties) for tensor products of homogeneous vector bundles and Borel–Weil–Bott computations. We describe how to locally run them here below.


## Getting started

### 1. Check Python

You need **Python 3.10 or newer**. If you do not have Python, install a current stable release from [python.org](https://www.python.org/downloads/).

Open a terminal (Terminal on macOS or Linux, or PowerShell on Windows) and type:

```text
python --version
```

On some computers the command is `python3` or `py`. If necessary, try `python3 --version` or `py --version`, and use whichever reports Python 3.10 or newer in place of `python` in all commands below.

### 2. Download the two repositories

Download [the-last-CY-pair](https://github.com/marcorampazzo/the-last-CY-pair) and [homogeneous-varieties](https://github.com/marcorampazzo/homogeneous-varieties) using **Code → Download ZIP** on each GitHub page. Extract both ZIP files and place the resulting folders inside the same folder. Rename them to remove any suffix such as `-main`, so that the layout is:

```text
some-folder/
├── homogeneous-varieties/
│   ├── grassmannians.py
│   └── ...
└── the-last-CY-pair/
    ├── common.py
    ├── hodge_numbers.py
    ├── pretilting.py
    └── open_cohomology.py
```

Keep the complete contents of both repositories. The scripts find the library automatically in this arrangement. Python and these two folders are all you need: there are no additional packages to install, and no virtual environment, conda, SageMath, or NumPy is required.

If you already use Git, you can instead run these two commands from the same folder:

```text
git clone https://github.com/marcorampazzo/homogeneous-varieties.git
git clone https://github.com/marcorampazzo/the-last-CY-pair.git
```

### 3. Run a first computation

In your terminal, enter the `the-last-CY-pair` folder. For example, replace the path below with its full location on your computer, keeping the quotation marks:

```text
cd "path/to/the-last-CY-pair"
```

For example, if you put both repositories inside a folder called `work` in your Downloads folder, use:

**macOS / Linux**

```text
cd "$HOME/Downloads/work/the-last-CY-pair"
```

**Windows PowerShell**

```text
cd "$HOME\Downloads\work\the-last-CY-pair"
```

`$HOME` automatically refers to your own user folder; you do not need to replace it with your name.

Then run:

```text
python pretilting.py
```

You should see:

```text
PASS: 225 ordered pairs; 654 finite cases; all remaining i >= j >= 0 by dominance.
Largest pair-specific tail threshold: i=6.
```

This verifies the higher-cohomology vanishings used in Lemma 6.8. The script checks 654 cases for the 225 ordered pairs of window bundles; the dominance argument in Appendix A.4 covers all remaining twists.

## Reproducing the paper

Run each command from the `the-last-CY-pair` folder. The scripts are independent and may be run in any order.

| Computation | Where it appears in the paper | Command |
| --- | --- | --- |
| Hodge diamonds of $X$ and $\widetilde Y$ | Section 2, Proposition 2.6, Appendix A.1 | `python hodge_numbers.py` |
| Higher-cohomology vanishings for the window bundles | Lemma 6.8, Appendix A.4 | `python pretilting.py` |
| Cohomology degrees of $\mathcal O(ah+bH)$ on $U/\Gamma$ | Corollary 6.18, Appendix A.5 | `python open_cohomology.py` |

### Hodge numbers

`hodge_numbers.py` prints both Hodge diamonds, with these values:

| Variety | $h^{1,1}$ | $h^{1,2}$ | Topological Euler characteristic |
| --- | --- | --- | --- |
| $X$ | 1 | 52 | −102 |
| $\widetilde Y$ (printed as `Ytilde`) | 2 | 53 | −102 |

The calculation assumes the smooth regular zero loci described in the paper. For $\widetilde Y$, it also uses the geometric fact that $h^{1,1}(\widetilde Y)\geq 2$: the blow-up supplies an exceptional divisor class in addition to a class pulled back from $Y$. As explained in Appendix A.1, this resolves the ambiguity left by the cohomological bounds. The script prints the bounds before applying this input. The Hodge numbers of $Y$ then follow from the blow-up formula, as in Proposition 2.6.

### Cohomology on the open quotient

`open_cohomology.py` prints a table with rows $-4\leq a\leq4$ and columns $-6\leq b\leq6$. Each entry lists the degrees in which $H^\bullet(U/\Gamma,\mathcal O(ah+bH))$ is nonzero. For example, `0,1,8` means nonzero cohomology in degrees 0, 1, and 8; the entries are degrees, not dimensions.

After the table, the script reports:

```text
PASS: all 13 ranges in the nine corollary items; no arbitrary truncation.
H0 is infinite dimensional for all (a,b); H1 is infinite dimensional exactly for a <= -2.
```

The infinite sums are handled by the finite reduction explained in Appendix A.5.


## Practical notes

The outputs above were checked with Python 3.12.14. On that run, the Hodge computation took about 8 seconds and each other command took less than a second; timings depend on your computer.

Results appear in the terminal. Each script checks its computed results and stops with an error if a check fails. Keep `common.py` beside the scripts: it contains shared setup and checking functions, and you do not need to run it.

If Python cannot find a script, check that your terminal is in `the-last-CY-pair`. If the library cannot be found, check the folder names and layout above; in particular, `grassmannians.py` must be directly inside `homogeneous-varieties`.
