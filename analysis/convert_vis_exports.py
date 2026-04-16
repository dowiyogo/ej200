from pathlib import Path
import subprocess
import sys


def main() -> int:
    src_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("vis_scan_png")
    if not src_dir.exists():
        print(f"ERROR: '{src_dir}' no existe.")
        return 1

    pdfs = sorted(src_dir.glob("trace_x_*.pdf"))
    if not pdfs:
        print(f"ERROR: no se encontraron PDFs en '{src_dir}'.")
        return 1

    for pdf in pdfs:
        png = pdf.with_suffix(".png")
        cmd = ["magick", "-density", "180", str(pdf), "-quality", "100", str(png)]
        print("Convirtiendo", pdf.name, "->", png.name)
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as exc:
            print(f"ERROR al convertir {pdf.name}: {exc}")
            return 1

    print("Listo.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
