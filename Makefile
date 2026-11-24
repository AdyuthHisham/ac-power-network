PY ?= python
.PHONY: figures
figures:
	cd analysis && $(PY) diagram.py && $(PY) e1_model.py && $(PY) e2_structure.py && $(PY) e3_lqg.py && $(PY) report_check.py
