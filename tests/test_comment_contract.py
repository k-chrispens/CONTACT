# - Locate authored Python source and test files.
from pathlib import Path

# - Verify short comments and bounded code sections.
def test_bullet_contract():
    root = Path(__file__).parents[1]
    files = [item for directory in ("src", "tests", "examples", "benchmarks")
             for item in (root / directory).rglob("*.py")]
    for item in files:
        count = 0
        assert item.read_text().splitlines()[0].startswith("# - "), item
        # - Count all executable lines between bullet comments.
        for number, line in enumerate(item.read_text().splitlines(), 1):
            text = line.lstrip()
            if text.startswith("# - "):
                assert len(text[4:].split()) <= 8, (item, number, "long comment")
                assert count <= 10, (item, number, "long code block", count)
                count = 0
            elif text and not text.startswith("#"):
                count += 1
        assert count <= 10, (item, "long final block", count)
