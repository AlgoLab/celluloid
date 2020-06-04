import unittest
import Translator
from tatsu.exceptions import FailedParse


class MyTestCase(unittest.TestCase):
    def test_EmptySASC(self):
        self.assertEqual(Translator.parse_file("testEmpty.txt", "SASC"), [])
    def test_EmptySCITE(self):
        self.assertEqual(Translator.parse_file("testEmpty.txt", "SCITE"), [])
    def test_EmptySPHYR(self):
        with self.assertRaises(FailedParse):
            Translator.parse_file("testEmpty.txt", "SPHYR")
    def test_SASC2SASC(self):
        ast1 = [['0', '1', '2', '0'], ['2', '0', '0', '0'], ['0', '1', '0', '0'], ['0', '1', '2', '0']]
        ast2 = [['0', '1', '2', '0'], ['2', '0', '0', '0'], ['0', '1', '0', '0'], ['0', '1', '2', '0']]
        ast3 = Translator.translate(ast1, "SASC", "SASC")
        self.assertEqual(ast3, ast2)
    def test_SASC2SPHYR(self):
        ast1 = [['0', '1', '2', '0'], ['2', '0', '0', '0'], ['0', '1', '0', '0'], ['0', '1', '2', '0']]
        ast2 = [['0', '1', '-1', '0'], ['-1', '0', '0', '0'], ['0', '1', '0', '0'], ['0', '1', '-1', '0']]
        ast3 = Translator.translate(ast1, "SASC", "SPHYR")
        self.assertEqual(ast3, ast2)
    def test_SASC2SCITE(self):
        ast1 = [['0', '1', '2', '0'], ['2', '0', '0', '0'], ['0', '1', '0', '0'], ['0', '1', '2', '0']]
        ast2 = [['0', '3', '0', '0'], ['1', '0', '1', '1'], ['2', '0', '0', '2'], ['0', '0', '0', '0']]
        ast3 = Translator.translate(ast1, "SASC", "SPHYR")
        self.assertEqual(ast3, ast2)


if __name__ == '__main__':
    unittest.main()
