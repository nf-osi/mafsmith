SRC      := manuscript.md
BIB      := references.bib
OUT_DIR  := build

PANDOC_FLAGS := --citeproc --bibliography=$(BIB) --standalone

.PHONY: all html pdf docx clean

all: html pdf docx

html: $(OUT_DIR)/manuscript.html
pdf:  $(OUT_DIR)/manuscript.pdf
docx: $(OUT_DIR)/manuscript.docx

$(OUT_DIR):
	mkdir -p $@

$(OUT_DIR)/manuscript.html: $(SRC) $(BIB) | $(OUT_DIR)
	pandoc $(PANDOC_FLAGS) -o $@ $(SRC)

$(OUT_DIR)/manuscript.pdf: $(SRC) $(BIB) | $(OUT_DIR)
	pandoc $(PANDOC_FLAGS) --pdf-engine=xelatex -o $@ $(SRC)

$(OUT_DIR)/manuscript.docx: $(SRC) $(BIB) | $(OUT_DIR)
	pandoc $(PANDOC_FLAGS) -o $@ $(SRC)

clean:
	rm -rf $(OUT_DIR)
