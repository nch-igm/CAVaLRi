from bs4 import BeautifulSoup

def parse_prettest_probability(table):

    # Drill down to second row of probability table
    trows = table.find_all('tr')
    content_row = trows[1]
    tcells = content_row.find_all('td')
    content_cell = tcells[0]
    prettest_prob = content_cell.text
    prettest_prob = prettest_prob.split('/')
    return float(prettest_prob[0])/float(prettest_prob[1])


def get_pheno_data(article):

    # Get pre-test probability
    probability_tables = article.find_all('table')
    for table in probability_tables:
        if "redTable" in table['class']:
            pretest_probability = parse_prettest_probability(table)


    # Get HPO-level data
    rects = article.find_all('rect')
    polygons = article.find_all('polygon')

    # Add polygon to rect list
    for p in polygons:
        rects.append(p)

    texts = article.find_all('text')
    hpo_ids = []
    for rect in rects:
        rect_info = rect['onmouseover']
        rect_info = rect_info[rect_info.find("'")+1:]

        # Check for an HPO ID listing inversion
        if rect_info.find('D<') == 0:
            start_pos = rect_info.find('<', 4)
        else:
            start_pos = 0

        if rect_info.find('[', start_pos) != -1:
            hpo_id = rect_info[rect_info.find('[', start_pos)+1:rect_info.find(']', start_pos)]
            lr = rect_info[rect_info.find('LR: ')+len('LR: '):rect_info.find("'", rect_info.find('LR: '))]
            hpo_ids.append({'hpoId': hpo_id, 'hpoLR': float(lr)})
    
    return hpo_ids


def get_variant_data(article):

    # Get the variant data
    tables = article.find_all('table')
    for table in tables:
        if "minimalistBlack" in table['class']:
            
            # Get each variant row
            variant_rows = table.findChildren('tr', recursive=False)

            for variant_row in variant_rows:    

                # Get gene and mode of inheritance
                if variant_row.text.find('Genotype score') != -1:

                    # Parse out gene and MOI
                    text = variant_row.find_all('td')[1].text
                    gene = text[:text.find(':')]
                    
                    # Get MOI
                    if text.find('Mode of inheritance:') != -1:
                        moi_pos = text.find('Mode of inheritance:') + len(' Mode of inheritance:')
                        moi = text[moi_pos: text.find('.', moi_pos)]
                    else:
                        moi = ''

                    # Get disease frequency
                    if text.find('disease=') != -1:
                        disease_freq_pos = text.find('disease=') + len('disease=')
                        disease_freq = float(text[disease_freq_pos: text.find('. ', disease_freq_pos)])
                    else:
                        disease_freq = 0.001

                    # Get background frequency
                    if text.find('background=') != -1:
                        background_freq_pos = text.find('background=') + len('background=')
                        background_freq = float(text[background_freq_pos: text.find('. ', background_freq_pos)])
                    else:
                        background_freq = 0.001

                    # Get geneLR
                    scraped_geneLR_pos = text.find('(LR)=') + len('(LR)=')
                    if text[-1] == '.':
                        scraped_geneLR = float(text[scraped_geneLR_pos:-1])
                    else:
                        scraped_geneLR = float(text[scraped_geneLR_pos:])

                    # Initialize omim data and add gene data
                    gene_data = {
                        'gene': gene,
                        'moi': moi,
                        'disease_freq': disease_freq,
                        'background_freq': background_freq,
                        'variants': []
                    }
                    
                
            for variant_row in variant_rows:

                # Only grab variant rows
                if variant_row.text.find('Genotype score') == -1:

                    # Get variant cells
                    variant_cells = variant_row.findChildren('td', recursive=False)

                    # Initialize counter
                    i = 0

                    for vc in variant_cells:
                        
                        # Get clinvar annotation
                        if i == 4:
                            clinvar = vc.text
                            i = 5

                        # Get genotype
                        if i == 3:
                            genotype = vc.text
                            i = 4

                        # Get population frequency
                        if i == 2:
                            popFreq = vc.text
                            i = 3

                        # Get Exomiser score
                        if i == 1:
                            pathScore = vc.text
                            i = 2
                        
                        # Get chromosome
                        if vc.text.find('chr') != -1:
                            pos = vc.text
                            i = 1

                    # Assemble data into dictionary
                    gene_data['variants'].append({
                        'pos': pos,
                        'pathScore': float(pathScore),
                        'popFreq': float(popFreq[: -1]),
                        'genotype': genotype,
                        'clinvar': clinvar
                        })

    gene_data.update({
        'calculated_geneLR': 0,#calculate_genotype_lr(gene_data),
        'scraped_geneLR': scraped_geneLR
    })

    return gene_data


def parse_lirical_html(case):

    with open(case.lirical_html_output) as f:
        soup = BeautifulSoup(f, 'html.parser')

    articles = soup.find_all('article')

    data = {
        'subjectId': case.case_id,
        'diseases': []
            }

    for article in articles:
        links = article.find_all('a')
        for link in links:
            try:
                if link['name'].find('diagnosis') != -1:

                    # Get disease name
                    article_name = link['name']
                    header = article.find_all('h3')
                    header = header[0].text
    
                    # Get OMIM ID
                    omim_id = header[header.find('[')+len('[OMIM:'):header.find(']')]
                    omim_name = header[header.find(' ') + 1:header.find('[')-1]
                    rank = int(header[header.find('(') + 1:header.find(')')])

                    # Get pre-test probability
                    probability_tables = article.find_all('table')
                    for table in probability_tables:
                        if "redTable" in table['class']:
                            pretest_probability = parse_prettest_probability(table)

                    # Get phenotype data
                    pheno_data = get_pheno_data(article)

                    # Get genotype data
                    gene_data = get_variant_data(article)
                   
                    # Append pheno data
                    data['diseases'].append({
                        'omimId': omim_id,
                        'diseaseName': omim_name,
                        'pretest_probability': pretest_probability,
                        'pheno_data': pheno_data,
                        'gene_data': gene_data
                    })

            except:
                pass
    
    return data


def build_case_data(case):
    return parse_lirical_html(case)
    