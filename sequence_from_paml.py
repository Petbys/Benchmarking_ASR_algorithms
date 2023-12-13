with open('/proj/naiss2023-6-162/applied_bioinformatics/aa_subset_Paml/rst', 'r') as file:
    text = file.read()

start_phrase = 'node #69'
end_phrase = 'node #70'

start_index = text.find(start_phrase)
end_index = text.find(end_phrase)

if start_index != -1 and end_index != -1:
    extracted_text = text[start_index:end_index]
    # Remove spaces from the extracted text
    extracted_text = extracted_text.replace(' ', '')
    print(extracted_text)
    