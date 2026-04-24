import unittest, importlib
from tempfile import TemporaryDirectory
import io
from contextlib import redirect_stdout
from find_errors import VariantWindow, find_errors


class Test_VariantWindow(unittest.TestCase):
	def test_add_gt(self):
		window = VariantWindow()
		with self.assertRaises(RuntimeError):
			window.add_gt("1|1")

		window.add_gt("0|1")
		self.assertTrue(window.nr_ref == 1);
		self.assertTrue(window.nr_alt == 0);
		self.assertTrue(window.genotypes == ["0|1"])

		window.add_gt("1|0")
		self.assertTrue(window.nr_ref == 1);
		self.assertTrue(window.nr_alt == 1);
		self.assertTrue(window.genotypes == ["0|1", "1|0"])

	def test_advance(self):
		window = VariantWindow()
		window.add_gt("0|1")
		window.add_gt("1|0")
		window.add_gt("0|1")
		window.add_gt("0|1")
		window.add_gt("0|1")

		self.assertTrue(window.nr_ref == 4)
		self.assertTrue(window.nr_alt == 1)
		self.assertTrue(window.genotypes == ["0|1", "1|0", "0|1", "0|1", "0|1"])

		window.advance("0|1")
		self.assertTrue(window.nr_ref == 4)
		self.assertTrue(window.nr_alt == 1)
		self.assertTrue(window.genotypes == ["1|0", "0|1", "0|1", "0|1", "0|1"])

		window.advance("0|1")
		self.assertTrue(window.nr_ref == 5)
		self.assertTrue(window.nr_alt == 0)
		self.assertTrue(window.genotypes == ["0|1", "0|1", "0|1", "0|1", "0|1"])
	
		window.advance("1|0")
		self.assertTrue(window.nr_ref == 4)
		self.assertTrue(window.nr_alt == 1)
		self.assertTrue(window.genotypes == ["0|1", "0|1", "0|1", "0|1", "1|0"])


	def test_check_position(self):
		window = VariantWindow()
		window.add_gt("0|1")
		window.add_gt("1|0")
		window.add_gt("0|1")
		window.add_gt("0|1")
		window.add_gt("0|1")	

		expected = [False, True, False, False, False]
		for i, e in enumerate(expected):
			self.assertTrue(window.check_position(i) == e)

		window = VariantWindow()
		window.add_gt("0|1")
		window.add_gt("1|0")
		window.add_gt("0|1")
		window.add_gt("1|0")
		window.add_gt("0|1")	

		expected = [False, True, False, True, False]
		for i, e in enumerate(expected):
			self.assertTrue(window.check_position(i) == e)

		window = VariantWindow()
		window.add_gt("1|0")
		window.add_gt("1|0")
		window.add_gt("0|1")
		window.add_gt("0|1")	

		expected = [False, False, False, False]
		for i, e in enumerate(expected):
			self.assertTrue(window.check_position(i) == e)

		window = VariantWindow()
		with self.assertRaises(RuntimeError):
			window.check_position(0)


class Test_find_errors(unittest.TestCase):
	def test_find_errors1(self):
		lines = ["chr1\t15952\t.\tT\tA\t.\tPASS\tAF=0.208333;UK=62;MA=0\tGT:KC\t1|0:2",
			 "chr1\t16509\t.\tA\tG\t.\tPASS\tAF=0.416667;UK=62;MA=0\tGT:KC\t1|1:3",
			 "chr1\t16636\t.\tT\tTA\t.\tPASS\tAF=0.375077;UK=50;MA=0\tGT:KC\t1|0:1",
			 "chr1\t18262\t.\tT\tC\t.\tPASS\tAF=0.916667;UK=62;MA=0\tGT:KC\t1|0:3"
		]

		expected_lines = ["chr1\t16509\t.\tA\tG\t.\tPASS\tAF=0.416667;UK=62;MA=0\tGT:KC\t1|1:3"]
		computed_lines = [l for l in find_errors(1, lines)]
		self.assertTrue(expected_lines == computed_lines)

	def test_find_errors2(self):
		lines = ["chr1\t15952\t.\tT\tA\t.\tPASS\tAF=0.208333;UK=62;MA=0\tGT:KC\t1|0:2",
			 "chr1\t16636\t.\tT\tTA\t.\tPASS\tAF=0.375077;UK=50;MA=0\tGT:KC\t0|1:1",
			 "chr1\t18262\t.\tT\tC\t.\tPASS\tAF=0.916667;UK=62;MA=0\tGT:KC\t1|0:3",
			 "chr1\t18262\t.\tT\tC\t.\tPASS\tAF=0.916667;UK=62;MA=0\tGT:KC\t0|1:3",
			 "chr1\t18262\t.\tT\tC\t.\tPASS\tAF=0.916667;UK=62;MA=0\tGT:KC\t1|0:3"
		]

		expected_lines = ["chr1\t16636\t.\tT\tTA\t.\tPASS\tAF=0.375077;UK=50;MA=0\tGT:KC\t0|1:1"]
		computed_lines = [l for l in find_errors(1, lines)]
		self.assertTrue(expected_lines == computed_lines)

	def test_find_errors3(self):
		lines = ["chr1\t15952\t.\tT\tA\t.\tPASS\tAF=0.208333;UK=62;MA=0\tGT:KC\t1|0:2",
			 "chr1\t16636\t.\tT\tTA\t.\tPASS\tAF=0.375077;UK=50;MA=0\tGT:KC\t1|0:1",
			 "chr1\t18262\t.\tT\tC\t.\tPASS\tAF=0.916667;UK=62;MA=0\tGT:KC\t1|0:3",
			 "chr1\t18262\t.\tT\tC\t.\tPASS\tAF=0.916667;UK=62;MA=0\tGT:KC\t0|1:3",
			 "chr1\t18262\t.\tT\tC\t.\tPASS\tAF=0.916667;UK=62;MA=0\tGT:KC\t0|1:3",
			"chr1\t18262\t.\tT\tC\t.\tPASS\tAF=0.916667;UK=62;MA=0\tGT:KC\t0|1:3"
		]

		expected_lines = []
		computed_lines = [l for l in find_errors(1, lines)]
		self.assertTrue(expected_lines == computed_lines)


	def test_find_errors4(self):
		lines = ["chr1\t15952\t.\tT\tA\t.\tPASS\tAF=0.208333;UK=62;MA=0\tGT:KC\t1|0:2",
			 "chr1\t16636\t.\tT\tTA\t.\tPASS\tAF=0.375077;UK=50;MA=0\tGT:KC\t1|0:1",
			 "chr1\t18262\t.\tT\tC\t.\tPASS\tAF=0.916667;UK=62;MA=0\tGT:KC\t1|0:3",
			 "chr1\t18262\t.\tT\tC\t.\tPASS\tAF=0.916667;UK=62;MA=0\tGT:KC\t0|1:3",
			 "chr1\t18262\t.\tT\tC\t.\tPASS\tAF=0.916667;UK=62;MA=0\tGT:KC\t0|1:3",
			"chr1\t18262\t.\tT\tC\t.\tPASS\tAF=0.916667;UK=62;MA=0\tGT:KC\t0|1:3"
		]

		expected_lines = []
		computed_lines = [l for l in find_errors(1, lines)]
		self.assertTrue(expected_lines == computed_lines)






if __name__ == '__main__':
    unittest.main()
