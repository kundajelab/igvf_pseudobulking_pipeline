import defopt

from visualize_qc.visualize_qc import visualize_qc


def main() -> None:
    """Run the visualize-qc command line interface."""
    defopt.run(visualize_qc)
