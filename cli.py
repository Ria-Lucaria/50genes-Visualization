import argparse
import sys

from combined import SUPPORTED_FORMATS, render_combined_figure


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Render combined figure from CSV.")
    parser.add_argument("-i", "--input", required=True, help="Path to input CSV file.")
    parser.add_argument("-o", "--output-dir", default="output", help="Output directory.")
    parser.add_argument("-p", "--prefix", default="Fig6", help="Output filename prefix.")
    parser.add_argument("-f", "--format", default="svg", choices=sorted(SUPPORTED_FORMATS), help="Output format.")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    try:
        output_file = render_combined_figure(
            csv_path=args.input,
            output_dir=args.output_dir,
            filename_prefix=args.prefix,
            output_format=args.format,
        )
        print(f"SUCCESS: {output_file}")
        return 0
    except (FileNotFoundError, ValueError, PermissionError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    except Exception as exc:  # pragma: no cover
        print(f"UNEXPECTED_ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
