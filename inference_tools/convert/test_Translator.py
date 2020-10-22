import unittest
import Translator
from tatsu.exceptions import FailedParse


class MyTestCase(unittest.TestCase):
    def setUp(self):
        self.formats  = ["SASC", "SPHYR", "SCITE"]
        self.ASTs = {
            "SASC" : [['0', '1', '2', '0'], ['2', '0', '0', '0'], ['0', '1', '0', '0'], ['0', '1', '2', '0']],
            "SPHYR" : [['0', '1', '-1', '0'], ['-1', '0', '0', '0'], ['0', '1', '0', '0'], ['0', '1', '-1', '0']],
            "SCITE" : [['0', '3', '0', '0'], ['1', '0', '1', '1'], ['3', '0', '0', '3'], ['0', '0', '0', '0']]
        }

    def test_readSASC(self):
        self.assertEqual(Translator.parse_string("", "SASC"), [])
        self.assertEqual(Translator.parse_string("0 1 2 0 \n 2 0 0 0 \n 0 1 0 0 \n 0 1 2 0", "SASC"),
                         [['0', '1', '2', '0'], ['2', '0', '0', '0'], ['0', '1', '0', '0'], ['0', '1', '2', '0']])
        with self.assertRaises(FailedParse):
            Translator.parse_string("0 1 2 0 \n 2 0 0 0 \n 0 1 0 0 \n 0 1 2 -1", "SASC")
        with self.assertRaises(Translator.NotAMatrix):
            Translator.parse_string("0 1 2 0 \n 2 0 0 0 \n 0 1 0 0 \n 0 1 2 ", "SASC")

    def test_readSCITE(self):
        self.assertEqual(Translator.parse_string("", "SCITE"), [])
        self.assertEqual(Translator.parse_string("0 3 0 0 \n 1 0 1 1 \n 3 0 0 3 \n 0 0 0 0 \n", "SCITE"),
                         [['0', '3', '0', '0'], ['1', '0', '1', '1'], ['3', '0', '0', '3'], ['0', '0', '0', '0']])
        with self.assertRaises(FailedParse):
            Translator.parse_string("0 3 0 0 \n 1 0 1 1 \n 3 0 0 3 \n 0 0 0 -1 \n", "SCITE")
        with self.assertRaises(Translator.NotAMatrix):
            Translator.parse_string("0 3 0 0 \n 1 0 1 1 \n 3 0 0 3 \n 0 0 0 \n", "SCITE")

    def test_readSPHYR(self):
        with self.assertRaises(FailedParse):
            Translator.parse_string("", "SPHYR")
        self.assertEqual(
            Translator.parse_string("4 #cells \n 4 #SNVs \n 0 1 -1 0 \n -1 0 0 0 \n 0 1 0 0 \n 0 1 -1 0", "SPHYR"),
            [['0', '1', '-1', '0'], ['-1', '0', '0', '0'], ['0', '1', '0', '0'], ['0', '1', '-1', '0']])
        with self.assertRaises(FailedParse):
            Translator.parse_string("0 1 2 0 \n 2 0 0 0 \n 0 1 0 0 \n 0 1 2 2", "SPHYR")
        with self.assertRaises(Translator.NotAMatrix):
            Translator.parse_string("4 #cells \n 4 #SNVs \n0 1 -1 0 \n -1 0 0 0 \n 0 1 0 0 \n 0 1", "SPHYR")


    def test_SASC(self):
        ast1 = [['0', '1', '2', '0'], ['2', '0', '0', '0'], ['0', '1', '0', '0'], ['0', '1', '2', '0']]
        for i in self.formats:
            ast3 = Translator.translate(ast1, "SASC", i)
            self.assertEqual(ast3, self.ASTs[i])

    def test_SPHYR(self):
        ast1 = [['0', '1', '-1', '0'], ['-1', '0', '0', '0'], ['0', '1', '0', '0'], ['0', '1', '-1', '0']]
        for i in self.formats:
            ast3 = Translator.translate(ast1, "SPHYR", i)
            self.assertEqual(ast3, self.ASTs[i])

    def test_SCITE(self):
        ast1 = [['0', '3', '0', '0'], ['1', '0', '2', '1'], ['3', '0', '0', '3'], ['0', '0', '0', '0']]
        for i in self.formats:
            ast3 = Translator.translate(ast1, "SCITE", i)
            self.assertEqual(ast3, self.ASTs[i])


if __name__ == '__main__':
    unittest.main()
