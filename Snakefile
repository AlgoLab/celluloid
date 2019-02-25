time = '/usr/bin/time'
_sasc_ = 'sasc/sasc'
_alpha_ = 0.25 # false negative rate
_beta_ = 0.00001 # false positive rate

methods = ['huang', 'random'] # initialization methods
dissims = ['conflict', 'matching'] # dissimilarity measures

# patterns
name_patt = '{name,[0-9A-z_]+}'
lab_patt = 'l{labels,(no|yes)}'
meth_patt = 'm{method,(' + '|'.join(methods) + ')}'
diss_patt = 'd{dissim,(' + '|'.join(dissims) + ')}'
pattern = name_patt + '.' + lab_patt + '.' + meth_patt + '.' + diss_patt + '.n{inits,[0-9]+}.k{k,[0-9]+}'
km_pattern = name_patt + '.k{k,[0-9]+}'

#
# master rule
#----------------------------------------------------------------------
rule master :
	input :
		# clustering the instances (Section 3.1)
		expand('data/exp{x}/sim_{id}.lno.m{method}.d{dissim}.n1.k100.celluloid/out',
			x = [1,2,3],
			id = range(1,51),
			method = methods,
			dissim = dissims),

		expand('data/exp{x}/sim_{id}.k100.kmeans/out',
			x = [1,2,3],
			id = range(1,51)),
	
		# run sasc on some clustered instances (Section 3.2)
		expand('data/exp{x}/sim_{id}.lno.m{method}.dconflict.n1.k100.celluloid/sim_{id}.sasc.out',
			x = [1,2,3],
			id = range(1,51),
			method = methods),

		expand('data/exp{x}/sim_{id}.k100.kmeans/sim_{id}.sasc.out',
			x = [1,2,3],
			id = range(1,51)),

		# run sasc on (clustered) mgh64 dataset (Section 3.3)
		'data/real/MGH64.lyes.mhuang.dconflict.n1.k100.celluloid/MGH64.sasc.out'

#
# run SASC on dataset (matrix, mutation labels) clustered /w celluloid
#----------------------------------------------------------------------
def labels_file(wildcards) :

	prefix = '{path}/{name}.l{labels}.m{method}.d{dissim}.n{inits}.k{k}.celluloid/'.format(**wildcards)
	suffix = 'LABELS_clustered.txt'
	if wildcards.labels == 'yes' :
		suffix = wildcards.name + '_snv-names_clustered.txt'

	return prefix + suffix

# run sasc
rule run_sasc_celluloid :
	output : '{path}/' + pattern + '.celluloid/{name}.sasc.out'

	input :
		prgm = _sasc_,
		badge = '{path}/' + pattern + '.celluloid/out'

	params :
		matrix = '{path}/' + pattern + '.celluloid/{name}_scs_clustered.txt',
		mutations = labels_file,
		a = _alpha_,
		b = _beta_

	log :
		log = '{path}/' + pattern + '.celluloid/{name}.sasc.log',
		time = '{path}/' + pattern + '.celluloid/{name}.sasc.time'

	run :
		with open(params.matrix) as mat :
			m = len(mat.readline().split())
			n = sum(1 for line in mat) + 1
		shell('''

   {time} -v -o {log.time} \
      {input.prgm} -i {params.matrix} -e {params.mutations} \
         -n {n} -m {m} \
         -a {params.a} -b {params.b} -k 0 -x -l \
            > {output} 2> {log.log} ''')

#
# run SASC on dataset (matrix, mutation labels) clustered /w k-means
#----------------------------------------------------------------------
rule run_sasc_kmeans :
	output : '{path}/' + km_pattern + '.kmeans/{name}.sasc.out'

	input :
		prgm = _sasc_,
		badge = '{path}/' + km_pattern + '.kmeans/out'

	params :
		matrix = '{path}/' + km_pattern + '.kmeans/{name}_scs_kmean.matrix',
		mutations = '{path}/' + km_pattern + '.kmeans/{name}_scs_kmean.mutations',
		a = _alpha_,
		b = _beta_,

	log :
		log = '{path}/' + km_pattern + '.kmeans/{name}.sasc.log',
		time = '{path}/' + km_pattern + '.kmeans/{name}.sasc.time'

	run :
		with open(params.matrix) as mat :
			m = len(mat.readline().split())
			n = sum(1 for line in mat) + 1
		shell('''

   {time} -v -o {log.time} \
      {input.prgm} -i {params.matrix} -e {params.mutations} \
         -n {n} -m {m} \
         -a {params.a} -b {params.b} -k 0 -x -l \
            > {output} 2> {log.log} ''')

#
# cluster a dataset (matrix, mutation labels) with celluloid
#----------------------------------------------------------------------
rule celluloid :
	output : '{path}/' + pattern + '.celluloid/out'

	input :
		prgm = 'celluloid.py',
		matrix = '{path}/{name}_scs.txt'

	params :
		outdir = '{path}/' + pattern + '.celluloid',
		labels = lambda wildcards :
			'--labels {path}/{name}_snv-names.txt'.format(**wildcards) if wildcards.labels == 'yes' else ''

	log :
		log = '{path}/' + pattern + '.celluloid/log',
		time = '{path}/' + pattern + '.celluloid/time'

	shell : '''

   {time} -v -o {log.time} \
      python3 {input.prgm} -i {input.matrix} -o {params.outdir} -v \
         -m {wildcards.method} -d {wildcards.dissim} \
         -n {wildcards.inits} -k {wildcards.k} {params.labels} \
            > {log.log} 2>&1
   touch {output} '''

#
# cluster a matrix with k-means
#----------------------------------------------------------------------
rule kmeans :
	output : '{path}/' + km_pattern + '.kmeans/out'

	input :
		prgm = 'kmeans.py',
		matrix = '{path}/{name}_scs.txt'

	params :
		outdir = '{path}/' + km_pattern + '.kmeans'

	log :
		log = '{path}/' + km_pattern + '.kmeans/log',
		time = '{path}/' + km_pattern + '.kmeans/time'

	shell : '''

   {time} -v -o {log.time} \
      python3 {input.prgm} -f {input.matrix} -o {params.outdir} \
         -k {wildcards.k} -c mutations \
            > {log.log} 2>&1
   touch {output} '''
