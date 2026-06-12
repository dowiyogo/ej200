#!/usr/bin/env python3
"""Build a labeled contact sheet from rendered slide PNG files."""
import argparse
import math
import pathlib

from PIL import Image, ImageDraw


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("rendered", type=pathlib.Path)
    parser.add_argument("output", type=pathlib.Path)
    parser.add_argument("--columns", type=int, default=5)
    args = parser.parse_args()

    slides = sorted(args.rendered.glob("slide-*.png"))
    if not slides:
        raise SystemExit(f"no rendered slides in {args.rendered}")
    thumb_w, thumb_h, label_h = 320, 180, 24
    rows = math.ceil(len(slides) / args.columns)
    sheet = Image.new("RGB", (args.columns * thumb_w, rows * (thumb_h + label_h)), "white")
    draw = ImageDraw.Draw(sheet)
    for index, slide in enumerate(slides):
        image = Image.open(slide).convert("RGB")
        image.thumbnail((thumb_w, thumb_h))
        x = (index % args.columns) * thumb_w
        y = (index // args.columns) * (thumb_h + label_h)
        sheet.paste(image, (x + (thumb_w - image.width) // 2, y))
        draw.text((x + 6, y + thumb_h + 3), f"slide {index + 1}", fill="black")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    sheet.save(args.output)


if __name__ == "__main__":
    main()
