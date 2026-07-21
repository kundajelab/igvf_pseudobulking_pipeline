import defopt

from visualize_qc.visualize_qc import visualize_qc


def main() -> None:
    defopt.run(visualize_qc)
