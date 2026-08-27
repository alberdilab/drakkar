import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "drakkar" / "workflow" / "scripts" / "validate_card_database.py"
SPEC = importlib.util.spec_from_file_location("validate_card_database", SCRIPT)
validator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(validator)


class ValidateCardDatabaseTests(unittest.TestCase):
    def test_reads_embedded_card_version(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            card_json = Path(tmpdir) / "card.json"
            card_json.write_text(json.dumps({"1": {}, "_version": "4.0.2"}), encoding="utf-8")
            self.assertEqual(validator.card_version(card_json), "4.0.2")

    def test_missing_version_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            card_json = Path(tmpdir) / "card.json"
            card_json.write_text("{}", encoding="utf-8")
            with self.assertRaises(ValueError):
                validator.card_version(card_json)


if __name__ == "__main__":
    unittest.main()
