#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <utility>
#include <algorithm>
#include <string>
#include <random>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <omp.h>
#include <fstream>
#include <map>

using namespace std;

static size_t signatureWidth; // Signature size (in bits)
static size_t signatureSize;  // Signature size (in uint64_t)
static size_t kmerLength;     // Kmer length
static size_t kmean_k;		  // Kmean k-value
static float density;         // % of sequence set as bits
static bool fastaOutput;      // Output fasta or csv

vector<pair<string, string>> loadFasta(const char *path)
{
	vector<pair<string, string>> sequences;

	FILE *fp = fopen(path, "r");
	if (!fp) {
		fprintf(stderr, "Failed to load %s\n", path);
		exit(1);
	}
	for (;;) {
		char seqNameBuf[8192];
		if (fscanf(fp, " >%[^\n]\n", seqNameBuf) < 1) break;
		string sequenceBuf;

		for (;;) {
			int c = fgetc(fp);
			if (c == EOF || c == '>') {
				ungetc(c, fp);
				break;
			}
			if (isalpha(c)) {
				sequenceBuf.push_back(c);
			}
		}
		sequences.push_back(make_pair(string(seqNameBuf), sequenceBuf));
	}
	fclose(fp);

	return sequences;
}

vector<uint64_t> readSignatures(const string file)
{
	ifstream rf(file, ios::out | ios::binary);
	// get length of file:
	rf.seekg(0, rf.end);
	int length = rf.tellg() / sizeof(uint64_t);
	rf.seekg(0, rf.beg);

	vector<uint64_t> sigs(length);
	size_t i = 0;
	while (rf) {
		rf.read((char *)&sigs[i], sizeof(uint64_t));
		i++;
	}
	rf.close();

	return sigs;
}

void generateSignature(uint64_t *output, const pair<string, string> &fasta)
{
	// Generate a signature from the kmers contained within

	string fastaSequence = fasta.second;
	// If the sequence is shorter than the kmer length, pad it with Xs
	while (fastaSequence.size() < kmerLength) {
		fastaSequence.push_back('X');
	}

	ranlux24_base rng;
	uniform_int_distribution<int> dist(-64 * signatureSize, signatureSize * 64 - 1);
	vector<int> unflattenedSignature(signatureSize * 64);
	int setBits = density * signatureSize * 64;
	//fprintf(stderr, "%s\n", fastaSequence.c_str());

	for (size_t i = 0; i < fastaSequence.size() - kmerLength + 1; i++) {
		seed_seq rngSeed(begin(fastaSequence) + i, begin(fastaSequence) + i + kmerLength);
		rng.seed(rngSeed);
		string kmer(begin(fastaSequence) + i, begin(fastaSequence) + i + kmerLength);
		//fprintf(stderr, "- %s\n", kmer.c_str());

		for (int j = 0; j < setBits; j++) {
			int bitPos = dist(rng);
			if (bitPos >= 0) {
				unflattenedSignature[bitPos] += 1;
			}
			else {
				unflattenedSignature[bitPos + 64 * signatureSize] -= 1;
			}
		}
	}
	fill(output, output + signatureSize, 0);
	for (size_t i = 0; i < signatureSize * 64; i++) {
		if (unflattenedSignature[i] > 0) {
			output[i / 64] |= (uint64_t)1 << (i % 64);
		}
	}
}

vector<uint64_t> convertFastaToSignatures(const vector<pair<string, string>> &fasta)
{
	vector<uint64_t> output;
	// Allocate space for the strings

	output.resize(fasta.size() * signatureSize);

#pragma omp parallel for schedule(dynamic)
	for (size_t i = 0; i < fasta.size(); i++) {
		generateSignature(&output[signatureSize * i], fasta[i]);
	}

	return output;
}

vector<uint64_t> filterSignatures(vector<uint64_t> signatures)
{
	map<vector<uint64_t>, size_t> filteredMap;
	for (size_t i = 0; i < signatures.size() / signatureSize; i++) {
		vector<uint64_t> sig(signatureSize);
		memcpy(&sig[0], &signatures[signatureSize * i], signatureSize * sizeof(uint64_t));
		filteredMap[sig]++;
	}

	vector<uint64_t> output(filteredMap.size()*signatureSize);
	size_t i = 0;
	for (map<vector<uint64_t>, size_t>::iterator it = filteredMap.begin(); it != filteredMap.end(); ++it) {
		memcpy(&output[signatureSize * i], &it->first[0], signatureSize * sizeof(uint64_t));
		i++;
	}
	return output;
}

void outputClusters(const vector<size_t> &clusters)
{
	for (size_t sig = 0; sig < clusters.size(); sig++)
	{
		printf("%llu,%llu\n", static_cast<unsigned long long>(sig), static_cast<unsigned long long>(clusters[sig]));
	}
}

void outputClusters(FILE * pFile, const vector<size_t> &clusters)
{
	for (size_t sig = 0; sig < clusters.size(); sig++)
	{
		fprintf(pFile, "%llu,%llu\n", static_cast<unsigned long long>(sig), static_cast<unsigned long long>(clusters[sig]));
	}
}

void outputFastaClusters(const vector<size_t> &clusters, const vector<pair<string, string>> &fasta)
{
	fprintf(stderr, "Writing out %zu records\n", clusters.size());
	for (size_t sig = 0; sig < clusters.size(); sig++)
	{
		printf(">%llu\n%s\n", static_cast<unsigned long long>(clusters[sig]), fasta[sig].second.c_str());
	}
}
/*
vector<size_t> clusterSignatures(const vector<uint64_t> &sigs)
{
auto rng = ranlux24_base();

auto dist = uniform_int_distribution<size_t>(0, clusterCount - 1);
size_t sigCount = sigs.size() / signatureSize;
vector<size_t> clusters(sigCount);

for (size_t i = 0; i < sigCount; i++) {
clusters[i] = dist(rng);
}

return clusters;
}
*/

// Parameters
size_t ktree_order = 10;
size_t ktree_capacity = 1000000;

void dbgPrintSignature(const uint64_t *sig)
{
	//fprintf(stderr, "%p: ", sig);
	for (size_t i = 0; i < signatureSize * 64; i++) {
		if (sig[i / 64] & (1ull << (i % 64))) {
			fprintf(stderr, "1");
		}
		else {
			fprintf(stderr, "0");
		}
	}
	fprintf(stderr, "\n");
}

void dbgPrintMatrix(const uint64_t *matrix)
{
	size_t ktree_csig_height = (ktree_order + 63) / 64;
	for (size_t i = 0; i < signatureSize * 64; i++) {
		fprintf(stderr, "%03zu:", i);
		for (size_t j = 0; j < ktree_csig_height * 64; j++) {
			auto val = matrix[i * ktree_csig_height + (j / 64)];
			if (val & (1ull << (j % 64))) {
				fprintf(stderr, "1");
			}
			else {
				fprintf(stderr, "0");
			}
		}
		fprintf(stderr, "\n");
		if (i >= 5) {
			fprintf(stderr, "...............\n");
			break;
		}
	}
}

template<class RNG>
vector<uint64_t> createRandomSigs(RNG &&rng, const vector<uint64_t> &sigs)
{
	constexpr size_t clusterCount = 2;
	vector<uint64_t> clusterSigs(signatureSize * clusterCount);
	size_t signatureCount = sigs.size() / signatureSize;
	uniform_int_distribution<size_t> dist(0, signatureCount - 1);
	bool finished = false;

	unordered_set<string> uniqueSigs;
	for (size_t i = 0; i < signatureCount; i++) {
		size_t sig = dist(rng);
		string sigData(signatureSize * sizeof(uint64_t), ' ');
		memcpy(&sigData[0], &sigs[sig * signatureSize], signatureSize * sizeof(uint64_t));
		uniqueSigs.insert(sigData);
		if (uniqueSigs.size() >= clusterCount) {
			finished = true;
			break;
		}
	}

	size_t i = 0;
	for (const auto &sig : uniqueSigs) {
		memcpy(&clusterSigs[i * signatureSize], sig.data(), signatureSize * sizeof(uint64_t));
		i++;
	}

	if (!finished) {
		if (uniqueSigs.size() != 1) {
			fprintf(stderr, "This should not happen\n");
			exit(1);
		}
		for (size_t i = 0; i < signatureSize; i++) {
			clusterSigs.push_back(clusterSigs[i]);
		}
	}

	return clusterSigs;
}

vector<vector<size_t>> createClusterLists(const vector<size_t> &clusters)
{
	constexpr size_t clusterCount = 2;
	vector<vector<size_t>> clusterLists(clusterCount);
	for (size_t i = 0; i < clusters.size(); i++) {
		clusterLists[clusters[i]].push_back(i);
	}
	return clusterLists;
}

vector<vector<size_t>> createClusterLists(const vector<size_t> &clusters, size_t clusterCount)
{
	vector<vector<size_t>> clusterLists(clusterCount);
	for (size_t i = 0; i < clusters.size(); i++) {
		clusterLists[clusters[i]].push_back(i);
	}
	return clusterLists;
}


vector<uint64_t> createClusterSigs(const vector<vector<size_t>> &clusterLists, const vector<uint64_t> &sigs)
{
	constexpr size_t clusterCount = 2;
	vector<uint64_t> clusterSigs(signatureSize * clusterCount);
	//#pragma omp parallel
	{
		vector<int> unflattenedSignature(signatureWidth);
		//#pragma omp for
		for (size_t cluster = 0; cluster < clusterLists.size(); cluster++) {
			fill(begin(unflattenedSignature), end(unflattenedSignature), 0);

			for (size_t signature : clusterLists[cluster]) {
				const uint64_t *signatureData = &sigs[signatureSize * signature];
				for (size_t i = 0; i < signatureWidth; i++) {
					uint64_t signatureMask = (uint64_t)1 << (i % 64);
					if (signatureMask & signatureData[i / 64]) {
						unflattenedSignature[i] += 1;
					}
					else {
						unflattenedSignature[i] -= 1;
					}
				}
			}

			uint64_t *flattenedSignature = &clusterSigs[cluster * signatureSize];
			for (size_t i = 0; i < signatureWidth; i++) {
				if (unflattenedSignature[i] > 0) {
					flattenedSignature[i / 64] |= (uint64_t)1 << (i % 64);
				}
			}
		}
	}
	return clusterSigs;
}

vector<uint64_t> createClusterSigs(const vector<vector<size_t>> &clusterLists, const vector<uint64_t> &sigs, size_t clusterCount)
{
	vector<uint64_t> clusterSigs(signatureSize * clusterCount);
	//#pragma omp parallel
	{
		vector<int> unflattenedSignature(signatureWidth);
		//#pragma omp for
		for (size_t cluster = 0; cluster < clusterLists.size(); cluster++) {
			fill(begin(unflattenedSignature), end(unflattenedSignature), 0);

			for (size_t signature : clusterLists[cluster]) {
				const uint64_t *signatureData = &sigs[signatureSize * signature];
				for (size_t i = 0; i < signatureWidth; i++) {
					uint64_t signatureMask = (uint64_t)1 << (i % 64);
					if (signatureMask & signatureData[i / 64]) {
						unflattenedSignature[i] += 1;
					}
					else {
						unflattenedSignature[i] -= 1;
					}
				}
			}

			uint64_t *flattenedSignature = &clusterSigs[cluster * signatureSize];
			for (size_t i = 0; i < signatureWidth; i++) {
				if (unflattenedSignature[i] > 0) {
					flattenedSignature[i / 64] |= (uint64_t)1 << (i % 64);
				}
			}
		}
	}
	return clusterSigs;
}


void reclusterSignatures(vector<size_t> &clusters, const vector<uint64_t> &meanSigs, const vector<uint64_t> &sigs)
{
	set<size_t> allClusters;
	for (size_t sig = 0; sig < clusters.size(); sig++) {
		const uint64_t *sourceSignature = &sigs[sig * signatureSize];
		size_t minHdCluster = 0;
		size_t minHd = numeric_limits<size_t>::max();

		for (size_t cluster = 0; cluster < 2; cluster++) {
			const uint64_t *clusterSignature = &meanSigs[cluster * signatureSize];
			size_t hd = 0;
			for (size_t i = 0; i < signatureSize; i++) {
				hd += __builtin_popcountll(sourceSignature[i] ^ clusterSignature[i]);
			}
			if (hd < minHd) {
				minHd = hd;
				minHdCluster = cluster;
			}
		}
		clusters[sig] = minHdCluster;
		allClusters.insert(minHdCluster);
	}

	if (allClusters.size() == 1) {
		// We can't have everything in the same cluster.
		// If this did happen, just split them evenly
		for (size_t sig = 0; sig < clusters.size(); sig++) {
			clusters[sig] = sig % 2;
		}
	}
}

void reclusterSignatures(vector<size_t> &clusters, const vector<uint64_t> &meanSigs, const vector<uint64_t> &sigs, size_t clusterCount)
{
	set<size_t> allClusters;
	for (size_t sig = 0; sig < clusters.size(); sig++) {
		const uint64_t *sourceSignature = &sigs[sig * signatureSize];
		size_t minHdCluster = 0;
		size_t minHd = numeric_limits<size_t>::max();

		for (size_t cluster = 0; cluster < clusterCount; cluster++) {
			const uint64_t *clusterSignature = &meanSigs[cluster * signatureSize];
			size_t hd = 0;
			for (size_t i = 0; i < signatureSize; i++) {
				hd += __builtin_popcountll(sourceSignature[i] ^ clusterSignature[i]);
			}
			if (hd < minHd) {
				minHd = hd;
				minHdCluster = cluster;
			}
		}
		clusters[sig] = minHdCluster;
		allClusters.insert(minHdCluster);
	}

	if (allClusters.size() == 1) {
		// We can't have everything in the same cluster.
		// If this did happen, just split them evenly
		for (size_t sig = 0; sig < clusters.size(); sig++) {
			clusters[sig] = sig % clusterCount;
		}
	}
}

vector<size_t> compressClusterRMSD(vector<size_t> &clusters, vector<size_t> RMSDs)
{
	vector<size_t>new_rmsd;
	unordered_map<size_t, size_t> remap;
	vector<size_t> clusters2 = clusters;

	for (size_t &clus : clusters) {
		if (remap.count(clus)) {
			clus = remap[clus];
		}
		else {
			size_t newClus = remap.size();
			remap[clus] = newClus;

			size_t size = count(clusters2.begin(), clusters2.end(), clus);
			new_rmsd.push_back(RMSDs[clus]);

			clus = newClus;

		}
	}
	fprintf(stderr, "Output %zu clusters\n", remap.size());
	return new_rmsd;
}

// There are two kinds of ktree nodes- branch nodes and leaf nodes
// Both contain a signature matrix, plus their own signature
// (the root node signature does not matter and can be blank)
// Branch nodes then contain 'order' links to other nodes
// Leaf nodes do not.
// However, as leaf nodes may become branch nodes, we allocate
// the space anyway.
// As the space to be used is determined at runtime, we use
// parallel arrays, not structs

struct KTree {
	size_t root = numeric_limits<size_t>::max(); // # of root node
	vector<size_t> childCounts; // n entries, number of children
	vector<int> isBranchNode; // n entries, is this a branch node
							  //vector<size_t> childLinks; // n * o entries, links to children
	vector<vector<size_t>> childLinks; // n entries, lists of links to children
	vector<size_t> parentLinks; // n entries, links to parents
	vector<uint64_t> means; // n * signatureSize entries, node signatures
							//vector<uint64_t> matrices; // n * (o / 64) * signatureSize * 64 entries
	vector<vector<uint64_t>> sigList; // n entries,  list of signatures
	vector<omp_lock_t> locks; // n locks
	size_t order;
	size_t capacity = 0; // Set during construction, currently can't change
	size_t matrixHeight;
	size_t matrixSize;

	void reserve(size_t capacity) {
		// For safety, only call this at startup currently
		if (this->capacity != 0) {
			fprintf(stderr, "Reserve can only be called from 0 capacity\n");
			exit(1);
		}
		this->capacity = capacity;
		matrixHeight = (order + 63) / 64;
		matrixSize = matrixHeight * signatureSize * 64;

#pragma omp parallel
		{
#pragma omp single
		{
			childCounts.resize(capacity);
		}
#pragma omp single
		{
			isBranchNode.resize(capacity);
		}
		//#pragma omp single
		//		{
		//			childLinks.resize(capacity * order);
		//		}
#pragma omp single
		{
			childLinks.resize(capacity);
		}
#pragma omp single
		{
			parentLinks.resize(capacity);
		}
#pragma omp single
		{
			locks.resize(capacity);
		}
		//#pragma omp single
		//		{
		//			matrices.resize(capacity * matrixSize);
		//		}
#pragma omp single
		{
			sigList.resize(capacity);
		}
#pragma omp single
		{
			means.resize(capacity * signatureSize);
		}
		}
	}

	KTree(size_t order_, size_t capacity) : order{ order_ } {
		reserve(capacity);
	}

	size_t calcHD(const uint64_t *a, const uint64_t *b) const
	{
		size_t c = 0;
		for (size_t i = 0; i < signatureSize; i++) {
			c += __builtin_popcountll(a[i] ^ b[i]);
		}
		return c;
	}

	size_t calcRMSD(size_t node) {
		size_t sumSquareHD = 0;
		size_t children = childCounts[node];

		if (children == 0) {
			return 0;
		}

		for (size_t i = 0; i < children; i++) {
			const uint64_t *signatureData = &sigList[node][signatureSize * i];
			size_t HD = calcHD(&means[node * signatureSize], signatureData);
			sumSquareHD += HD * HD;
		}

		return sqrt(sumSquareHD / children);
	}

	size_t traverse(const uint64_t *signature) const
	{
		size_t node = root;
		while (isBranchNode[node]) {
			size_t lowestHD = numeric_limits<size_t>::max();
			size_t lowestHDchild = 0;

			for (size_t i = 0; i < childCounts[node]; i++) {
				//size_t child = childLinks[node * order + i];
				size_t child = childLinks[node][i];
				size_t hd = calcHD(&means[child * signatureSize], signature);
				if (hd < lowestHD) {
					lowestHD = hd;
					lowestHDchild = child;
				}
			}
			node = lowestHDchild;
		}
		return node;
	}

	void addSigToMatrix(uint64_t *matrix, size_t child, const uint64_t *sig) const
	{
		size_t childPos = child / 64;
		size_t childOff = child % 64;

		//fprintf(stderr, "Adding this signature:\n");
		//dbgPrintSignature(sig);
		//fprintf(stderr, "To this matrix:\n");
		//dbgPrintMatrix(matrix);

		for (size_t i = 0; i < signatureSize * 64; i++) {
			matrix[i * matrixHeight + childPos] |= ((sig[i / 64] >> (i % 64)) & 0x01) << childOff;
		}
		//fprintf(stderr, "Resulting in:\n");
		//dbgPrintMatrix(matrix);
	}
	void removeSigFromMatrix(uint64_t *matrix, size_t child) const
	{
		size_t childPos = child / 64;
		size_t childOff = child % 64;

		uint64_t mask = ~(1ull << childOff);

		//fprintf(stderr, "Removing the %zuth child from matrix\n", child);    
		for (size_t i = 0; i < signatureSize * 64; i++) {
			matrix[i * matrixHeight + childPos] &= mask;
		}
		//fprintf(stderr, "Resulting in:\n");
		//dbgPrintMatrix(matrix);
	}

	void addSigToSigList(size_t node, const uint64_t *sig) {
		for (size_t i = 0; i < signatureSize; i++) {
			sigList[node].push_back(sig[i]);
		}
	}

	//void recalculateSig(size_t node)
	//{
	//	size_t children = childCounts[node];
	//	uint64_t *matrix = &matrices[node * matrixSize];
	//	uint64_t *sig = &means[node * signatureSize];
	//	fill(sig, sig + signatureSize, 0ull);

	//	auto threshold = (children / 2) + 1;

	//	for (size_t i = 0; i < signatureSize * 64; i++) {
	//		size_t c = 0;
	//		for (size_t j = 0; j < matrixHeight; j++) {
	//			auto val = matrix[i * matrixHeight + j];
	//			c += __builtin_popcountll(val);
	//		}
	//		if (c >= threshold) {
	//			sig[i / 64] |= 1ull << (i % 64);
	//		}
	//	}
	//	//fprintf(stderr, "Mean sig:\n");
	//	//dbgPrintSignature(sig);
	//}


	void recalculateSig(size_t node)
	{
		uint64_t *meanSig = &means[node*signatureSize];
		fill(meanSig, meanSig + signatureSize, 0ull);
		vector<int> unflattenedSignature(signatureWidth);
		fill(begin(unflattenedSignature), end(unflattenedSignature), 0);

		for (size_t i = 0; i < childCounts[node]; i++) {
			const uint64_t *signatureData = &sigList[node][i*signatureSize];

			for (size_t i = 0; i < signatureWidth; i++) {
				uint64_t signatureMask = (uint64_t)1 << (i % 64);
				if (signatureMask & signatureData[i / 64]) {
					unflattenedSignature[i] += 1;
				}
				else {
					unflattenedSignature[i] -= 1;
				}
			}
		}

		// update node mean
		for (size_t i = 0; i < signatureWidth; i++) {
			if (unflattenedSignature[i] > 0) {
				meanSig[i / 64] |= (uint64_t)1 << (i % 64);
			}
		}

	}

	void recalculateUp(size_t node)
	{
		size_t limit = 10;
		//fprintf(stderr, "RecalculateUp %zu\n", node);
		while (node != root) {
			recalculateSig(node);

			// update parent matrix with new mean
			size_t parent = updateParentMatrix(node, &means[node * signatureSize]);
			omp_unset_lock(&locks[parent]);

			node = parentLinks[node];
			if (omp_test_lock(&locks[node])) {
				omp_unset_lock(&locks[node]);
			}
			else {
				break;
			}

			// Put a limit on how far we go up
			// At some point it stops mattering, plus this helps avoid inf loops
			// caused by cycles getting into the tree structure
			limit--;
			if (limit == 0) return;
			//fprintf(stderr, "-> %zu\n", node);
		}
	}

	size_t getNewNodeIdx(vector<size_t> &insertionList)
	{
		if (insertionList.empty()) {
			fprintf(stderr, "ERROR: ran out of insertion points\n");
			exit(1);
		}
		size_t idx = insertionList.back();
		insertionList.pop_back();

		// Initialise lock
		omp_init_lock(&locks[idx]);
		return idx;
	}

	// attempt to find parent
	size_t getParent(size_t node) {
		size_t parent = parentLinks[node];

		// Lock the parent
		omp_set_lock(&locks[parent]);

		size_t idx = numeric_limits<size_t>::max();
		for (size_t i = 0; i < childCounts[parent]; i++) {
			//if (childLinks[parent * order + i] == node) {
			if (childLinks[parent][i] == node) {
				idx = i;
				break;
			}
		}

		// get the wrong parent, try again
		size_t count = 0;
		while (idx == numeric_limits<size_t>::max()) {
			// unset old lock
			omp_unset_lock(&locks[parent]);

			// find and lock new parent
			parent = parentLinks[node];
			omp_set_lock(&locks[parent]);

			// find idx
			for (size_t i = 0; i < childCounts[parent]; i++) {
				//if (childLinks[parent * order + i] == node) {
				if (childLinks[parent][i] == node) {
					idx = i;
					break;
				}
			}
			count++;
			if (count > 10) { // try 10 times only
				omp_unset_lock(&locks[parent]);
				return -1;
			}
		}

		omp_unset_lock(&locks[parent]);
		return parent;
	}

	// lock parent while updating matrix, return parent node for unlocking later
	size_t updateParentMatrix(size_t node, uint64_t *meanSig) {
		size_t parent = parentLinks[node];

		// Lock the parent
		omp_set_lock(&locks[parent]);

		size_t idx = numeric_limits<size_t>::max();
		for (size_t i = 0; i < childCounts[parent]; i++) {
			//if (childLinks[parent * order + i] == node) {
			if (childLinks[parent][i] == node) {
				idx = i;
				break;
			}
		}

		// get the wrong parent, try again
		size_t count = 0;
		while (idx == numeric_limits<size_t>::max()) {
			// unset old lock
			omp_unset_lock(&locks[parent]);

			// find and lock new parent
			parent = parentLinks[node];
			omp_set_lock(&locks[parent]);

			// find idx
			for (size_t i = 0; i < childCounts[parent]; i++) {
				//if (childLinks[parent * order + i] == node) {
				if (childLinks[parent][i] == node) {
					idx = i;
					break;
				}
			}
			count++;
			// after 10 trials, skip update, return parent for unlocking
			if (count > 10) { 
				return parent;
			}
		}

		//removeSigFromMatrix(&matrices[parent * matrixSize], idx);
		//addSigToMatrix(&matrices[parent * matrixSize], idx, meanSig);

		memcpy(&sigList[parent][idx * signatureSize], meanSig, signatureSize * sizeof(uint64_t));

		return parent;
	}


	template<class RNG>
	void splitNode(RNG &&rng, size_t node, const uint64_t *sig, vector<size_t> &insertionList, size_t link)
	{
		//fprintf(stderr, "Splitting node %zu\n", node);
		// Add 'sig' to the current node, splitting it in the process
		//fprintf(stderr, "Adding signature:\n");
		//dbgPrintSignature(sig);
		size_t nodeSigs = childCounts[node] + 1;
		vector<uint64_t> sigs(nodeSigs * signatureSize);
		copy(sigList[node].begin(), sigList[node].end(), sigs.begin());
		memcpy(&sigs[childCounts[node] * signatureSize], sig, signatureSize * sizeof(uint64_t));



		//sigList[node].resize(nodeSigs * signatureSize);
		//memcpy(&sigList[node][childCounts[node] * signatureSize], sig, signatureSize * sizeof(uint64_t));
		//addSigToSigList(node, sig);

		/*for (int i = 0; i < childCounts[node]; i++) {
		uint64_t *currentSig = &sigs[i * signatureSize];
		uint64_t *matrix = &matrices[node * matrixSize];
		for (size_t j = 0; j < signatureSize * 64; j++) {
		currentSig[j / 64] |= ((matrix[j * matrixHeight + i / 64] >> (i % 64)) & 1) << (j % 64);
		}
		}*/


		/*fprintf(stderr, "Signatures converted for clustering:\n");
		for (size_t i = 0; i < nodeSigs; i++) {
		uint64_t *currentSig = &sigs[i * signatureSize];
		fprintf(stderr, "matrix\n");
		dbgPrintSignature(currentSig);
		}*/




		vector<uint64_t> meanSigs = createRandomSigs(rng, sigs);
		vector<size_t> clusters(nodeSigs);
		vector<vector<size_t>> clusterLists;
		for (int iteration = 0; iteration < 4; iteration++) {
			//fprintf(stderr, "Iteration %d\n", iteration);
			reclusterSignatures(clusters, meanSigs, sigs);
			clusterLists = createClusterLists(clusters);
			meanSigs = createClusterSigs(clusterLists, sigs);
		}



		//vector<uint64_t> meanSigs = createRandomSigs(rng, sigList[node]);
		//vector<size_t> clusters(nodeSigs);
		//vector<vector<size_t>> clusterLists;
		//for (int iteration = 0; iteration < 1; iteration++) {
		//	//fprintf(stderr, "Iteration %d\n", iteration);
		//	reclusterSignatures(clusters, meanSigs, sigList[node]);
		//	clusterLists = createClusterLists(clusters);
		//	meanSigs = createClusterSigs(clusterLists, sigList[node]);
		//}

		//// Display clusters (debugging purposes)
		//for (const auto &clusterList : clusterLists) {
		//fprintf(stderr, "Cluster:\n");
		//for (size_t seqIdx : clusterList) {
		//uint64_t *currentSig = &sigs[seqIdx * signatureSize];
		//dbgPrintSignature(currentSig);
		//}
		//}


		// Create the sibling node
		size_t sibling = getNewNodeIdx(insertionList);

		size_t newlyAddedIdx = childCounts[node];

		childCounts[sibling] = clusterLists[1].size();
		isBranchNode[sibling] = isBranchNode[node];
		//sigList[sibling].resize(childCounts[sibling] * signatureSize);
		{
			size_t siblingIdx = 0;
			for (size_t seqIdx : clusterLists[1]) {
				//if (seqIdx < newlyAddedIdx) {
				//	childLinks[sibling * order + siblingIdx] = childLinks[node * order + seqIdx];
				//}
				//else {
				//	childLinks[sibling * order + siblingIdx] = link;
				//}

				// If this is a branch node, relink the child to the new parent
				if (isBranchNode[sibling]) {
					if (seqIdx < newlyAddedIdx) {
						childLinks[sibling].push_back(childLinks[node][seqIdx]);
					}
					else {
						childLinks[sibling].push_back(link);
					}

					//parentLinks[childLinks[sibling * order + siblingIdx]] = sibling;
					parentLinks[childLinks[sibling][siblingIdx]] = sibling;
				}

				//addSigToMatrix(&matrices[sibling * matrixSize], siblingIdx, &sigs[seqIdx * signatureSize]);


				//memcpy(&sigList[sibling][siblingIdx * signatureSize], &sigs[seqIdx * signatureSize], signatureSize * sizeof(uint64_t));
				addSigToSigList(sibling, &sigs[seqIdx * signatureSize]);

				siblingIdx++;
			}
		}

		memcpy(&means[node * signatureSize], &meanSigs[0 * signatureSize], signatureSize * sizeof(uint64_t));
		memcpy(&means[sibling * signatureSize], &meanSigs[1 * signatureSize], signatureSize * sizeof(uint64_t));

		// Fill the current node with the other cluster of signatures
		childCounts[node] = clusterLists[0].size();
		{
			//childLinks[node].resize(childCounts[node]);
			//fill(&matrices[node * matrixSize], &matrices[node * matrixSize] + matrixSize, 0ull);
			sigList[node].clear();

			size_t nodeIdx = 0;
			for (size_t seqIdx : clusterLists[0]) {
				//if (seqIdx < newlyAddedIdx) {
				//	childLinks[node * order + nodeIdx] = childLinks[node * order + seqIdx];
				//}
				//else {
				//	childLinks[node * order + nodeIdx] = link;
				//}

				// If this is a branch node, relink the child to the new parent
				if (isBranchNode[node]) {
					if (seqIdx < newlyAddedIdx) {
						childLinks[node][nodeIdx] = childLinks[node][seqIdx];
					}
					else {
						childLinks[node][nodeIdx] = link;
					}

					//parentLinks[childLinks[node * order + nodeIdx]] = node;
					parentLinks[childLinks[node][nodeIdx]] = node;
				}
				//addSigToMatrix(&matrices[node * matrixSize], nodeIdx, &sigs[seqIdx * signatureSize]);

				//memcpy(&sigList[node][nodeIdx * signatureSize], &sigs[seqIdx * signatureSize], signatureSize * sizeof(uint64_t));
				addSigToSigList(node, &sigs[seqIdx * signatureSize]);
				nodeIdx++;
			}
			if (isBranchNode[node]) {
				childLinks[node].resize(childCounts[node]);
			}
		}
		//sigList[node].resize(childCounts[node] * signatureSize);

		// Is this the root level?
		if (node == root) {
			//fprintf(stderr, "Node being split is root node\n");

			// Create a new root node
			size_t newRoot;
			newRoot = getNewNodeIdx(insertionList);

			// Link this node and the sibling to it
			parentLinks[node] = newRoot;
			parentLinks[sibling] = newRoot;

			childCounts[newRoot] = 2;
			isBranchNode[newRoot] = 1;
			//childLinks[newRoot * order + 0] = node;
			//childLinks[newRoot * order + 1] = sibling;
			childLinks[newRoot].push_back(node);
			childLinks[newRoot].push_back(sibling);

			//addSigToMatrix(&matrices[newRoot * matrixSize], 0, &meanSigs[0 * signatureSize]);
			//addSigToMatrix(&matrices[newRoot * matrixSize], 1, &meanSigs[1 * signatureSize]);

			//sigList[newRoot].resize(2 * signatureSize);
			//memcpy(&sigList[newRoot][0 * signatureSize], &meanSigs[0 * signatureSize], signatureSize * sizeof(uint64_t));
			//memcpy(&sigList[newRoot][1 * signatureSize], &meanSigs[1 * signatureSize], signatureSize * sizeof(uint64_t));
			addSigToSigList(newRoot, &meanSigs[0 * signatureSize]);
			addSigToSigList(newRoot, &meanSigs[1 * signatureSize]);

			root = newRoot;
		}
		else {
			//// First, update the reference to this node in the parent with the new mean
			//size_t parent = parentLinks[node];

			//// Lock the parent
			//omp_set_lock(&locks[parent]);

			//size_t idx = numeric_limits<size_t>::max();
			//for (size_t i = 0; i < childCounts[parent]; i++) {
			//	if (childLinks[parent * order + i] == node) {
			//		idx = i;
			//		break;
			//	}
			//}
			//if (idx == numeric_limits<size_t>::max()) {
			//	fprintf(stderr, "Error: node %zu is not its parent's (%zu) child\n", node, parent);

			//	// Abort. Unlock the parent and get out of here
			//	omp_unset_lock(&locks[parent]);
			//	return;

			//	//exit(1);
			//}

			//removeSigFromMatrix(&matrices[parent * matrixSize], idx);
			//addSigToMatrix(&matrices[parent * matrixSize], idx, &meanSigs[0 * signatureSize]);

			// First, update the reference to this node in the parent with the new mean
			size_t parent = updateParentMatrix(node, &meanSigs[0 * signatureSize]);

			// Connect sibling node to parent
			parentLinks[sibling] = parent;

			// Now add a link in the parent node to the sibling node
			//size_t RMSD_parent = calcMatrixRMSD(parent, &matrices[parent*matrixSize], childCounts[parent]);
			size_t RMSD_parent = calcRMSD(parent);


			//if (RMSD_parent<RMSDthreshold){// || childCounts[parent] + 1 < 10){
			//if (RMSD_parent<threshold) {
			if (childCounts[parent] + 1 < order) {
				//addSigToMatrix(&matrices[parent * matrixSize], childCounts[parent], &meanSigs[1 * signatureSize]);


				//childLinks[parent * order + childCounts[parent]] = sibling;
				childLinks[parent].push_back(sibling);

				childCounts[parent]++;
				//sigList[parent].resize(childCounts[parent]);
				//memcpy(&sigList[parent][childCounts[parent] * signatureSize], &meanSigs[1 * signatureSize], signatureSize * sizeof(uint64_t));
				addSigToSigList(parent, &meanSigs[1 * signatureSize]);

				// Update signatures (may change?)
				recalculateUp(parent);
			}
			else {
				//fprintf(stderr, "parent %zu, node %zu, RMSD_parent %zu\n", parent, node, RMSD_parent);
				splitNode(rng, parent, &meanSigs[1 * signatureSize], insertionList, sibling);
			}
			// Unlock the parent
			omp_unset_lock(&locks[parent]);
		}

		//fprintf(stderr, "Split finished\n");
	}
	template<class RNG>
	void insert(RNG &&rng, const uint64_t *signature, vector<size_t> &insertionList)
	{
		// Warning: ALWAYS INSERT THE FIRST NODE SINGLE-THREADED
		// We don't have any protection from this because it would slow everything down to do so
		if (root == numeric_limits<size_t>::max()) {
			root = getNewNodeIdx(insertionList);
			childCounts[root] = 0;
			isBranchNode[root] = 0;
		}

		size_t insertionPoint = traverse(signature);

		//fprintf(stderr, "Inserting at %zu\n", insertionPoint);
		omp_set_lock(&locks[insertionPoint]);
		if (childCounts[insertionPoint] < order) {
			//addSigToMatrix(&matrices[insertionPoint * matrixSize], childCounts[insertionPoint], signature);
			addSigToSigList(insertionPoint, signature);
			childCounts[insertionPoint]++;
		}
		else {
			splitNode(rng, insertionPoint, signature, insertionList, 0);
		}
		omp_unset_lock(&locks[insertionPoint]);

		//fprintf(stderr, "Node %zu now has %zu leaves\n", insertionPoint, childCounts[insertionPoint]);
	}

	void destroyLocks(size_t node)
	{
		omp_destroy_lock(&locks[node]);
		if (isBranchNode[node]) {
			for (size_t i = 0; i < childCounts[node]; i++) {
				destroyLocks(childLinks[node][i]);
				//destroyLocks(childLinks[node * order + i]);
			}
		}
	}
	void destroyLocks()
	{
		destroyLocks(root);
	}

	void updateMeanSig(uint64_t *meanSig, const vector<uint64_t> &sigs, vector<size_t> sigIndices) {
		vector<int> unflattenedSignature(signatureWidth);
		fill(begin(unflattenedSignature), end(unflattenedSignature), 0);

		for (size_t signature : sigIndices) {
			const uint64_t *signatureData = &sigs[signatureSize * signature];

			for (size_t i = 0; i < signatureWidth; i++) {
				uint64_t signatureMask = (uint64_t)1 << (i % 64);
				if (signatureMask & signatureData[i / 64]) {
					unflattenedSignature[i] += 1;
				}
				else {
					unflattenedSignature[i] -= 1;
				}
			}
		}

		// update node mean
		fill(meanSig, meanSig + signatureSize, 0ull);
		for (size_t i = 0; i < signatureWidth; i++) {
			if (unflattenedSignature[i] > 0) {
				meanSig[i / 64] |= (uint64_t)1 << (i % 64);
			}
		}

	}

	size_t calcRMSD(uint64_t *meanSig, const vector<uint64_t> &sigs, vector<size_t> sigIndices) {
		size_t sumSquareHD = 0;
		size_t children = sigIndices.size();

		if (children == 0) {
			return 0;
		}

		for (size_t signature : sigIndices) {
			const uint64_t *signatureData = &sigs[signatureSize * signature];
			size_t HD = calcHD(meanSig, signatureData);
			sumSquareHD += HD * HD;
		}
		return sqrt(sumSquareHD / children);
	}

	// update node (and parent) mean, return list of RMSD
	vector<size_t> updateTree(vector<size_t> clusters, const vector<uint64_t> &sigs)
	{
		set<size_t> nonEmptyNodes(clusters.begin(), clusters.end());
		size_t maxClusterCount = *max_element(clusters.begin(), clusters.end()) + 1;
		vector<vector<size_t>> clusterLists = createClusterLists(clusters, maxClusterCount);
		vector<size_t> RMSDs(ktree_capacity);

		for (size_t node : nonEmptyNodes) {
			updateMeanSig(&means[node * signatureSize], sigs, clusterLists[node]);

			// update parent matrix
			size_t parent = updateParentMatrix(node, &means[node * signatureSize]);
			recalculateUp(parent);
			omp_unset_lock(&locks[parent]);

			// get distortion
			RMSDs[node] = calcRMSD(&means[node * signatureSize], sigs, clusterLists[node]);
		}
		return RMSDs;
	}

	template<class RNG>
	vector<size_t> kmeanCluster(RNG &&rng, const vector<uint64_t> &sigs, size_t clusterCount) {
		vector<uint64_t> meanSigs = createRandomSigs(rng, sigs);
		vector<size_t> clusters(sigs.size() / signatureSize);
		vector<vector<size_t>> clusterLists;
		for (size_t iteration = 0; iteration < 4; iteration++) {
			//printf(">\n");
			reclusterSignatures(clusters, meanSigs, sigs, clusterCount);
			clusterLists = createClusterLists(clusters, clusterCount);
			meanSigs = createClusterSigs(clusterLists, sigs, clusterCount);
		}
		return clusters;
	}

	// level 0 is bottom most nodes that holds the signatures
	size_t getNodesByLevel(size_t node, size_t level, set<size_t> &parents) {

		//fprintf(stderr, "Parent: %zu\n", node);
		for (size_t i = 0; i < childCounts[node]; i++) {
			size_t child = childLinks[node][i];
			if (isBranchNode[child]) {
				getNodesByLevel(child, level, parents);
			}
			else {
				size_t target = child;
				for (size_t j = 0; j < level; j++) {
					// root has no parent
					if (target == root) {
						return 0;
					}
					target = getParent(target);
				}
				//fprintf(stderr, "%zu\n", target);
				parents.insert(target);
			}
		}
	}

	template<class RNG>
	vector<size_t> clusterClusters(RNG &&rng, vector<size_t> inputClusters, const vector<uint64_t> &sigs) {
		// get meanSig of all leaf nodes
		set<size_t> nonEmptyNodes(inputClusters.begin(), inputClusters.end());
		size_t clusterCount = nonEmptyNodes.size();

		// use KTree cluster means as initial centroids
		vector<uint64_t> ktreeMeanSigs(clusterCount * signatureSize);

		size_t i = 0;
		for (size_t node : nonEmptyNodes) {
			memcpy(&ktreeMeanSigs[i * signatureSize], &means[node * signatureSize], signatureSize * sizeof(uint64_t));
			i++;
		}


		vector<uint64_t> meanSigs = kmeanCluster(rng, ktreeMeanSigs, kmean_k);

		vector<size_t> clusters(inputClusters.size());
		vector<vector<size_t>> clusterLists;
		//while (true) {
		int k_iteration = 20;
		for (size_t iteration = 0; iteration < k_iteration; iteration++) {
			//printf(">\n");
			//vector<size_t> clusters_temp = clusters;

			reclusterSignatures(clusters, meanSigs, sigs, kmean_k);
			clusterLists = createClusterLists(clusters, kmean_k);
			meanSigs = createClusterSigs(clusterLists, sigs, kmean_k);

			//iteration++;
			//string file_name = "silva-o" + to_string(ktree_order) +
			//	"-i" + to_string(iteration) + ".txt";
			//FILE * pFile = fopen(file_name.c_str(), "w");
			//outputClusters(pFile, clusters);

			//// reach convergence
			//if (clusters_temp == clusters) {
			//	fprintf(stderr, "Total iterations %d\n", iteration);
			//	break;
			//}
		}

		vector<size_t> RMSDs(ktree_capacity);
		for (size_t i = 0; i < clusterLists.size(); i++) {
			RMSDs[i] = calcRMSD(&meanSigs[i*signatureSize], sigs, clusterLists[i]);
		}

		vector<size_t> compressedRMSDs = compressClusterRMSD(clusters, RMSDs);

		// output kmean clusters
		string clus_file_name = "silva-o" + to_string(ktree_order) +
			"-i" + to_string(k_iteration) + ".txt";
		FILE * clusFile = fopen(clus_file_name.c_str(), "w");
		outputClusters(clusFile, clusters);

		//output kmean distortion
		string dist_file_name = "silva-o" + to_string(ktree_order) +
			"-i" + to_string(k_iteration) + "-distortion.txt";
		FILE * distFile = fopen(dist_file_name.c_str(), "w");
		for (size_t i = 0; i < compressedRMSDs.size(); i++) {
			fprintf(distFile, "%zu,%zu\n", i, compressedRMSDs[i]);
		}
		return clusters;
	}


	template<class RNG>
	vector<size_t> restructureTree(RNG &&rng, vector<size_t> inputClusters, const vector<uint64_t> &sigs) {
		// get meanSig of all leaf nodes
		set<size_t> nonEmptyNodes(inputClusters.begin(), inputClusters.end());
		size_t clusterCount = nonEmptyNodes.size();

		// use KTree cluster means as initial centroids
		vector<uint64_t> ktreeMeanSigs(clusterCount * signatureSize);

		size_t i = 0;
		for (size_t node : nonEmptyNodes) {
			memcpy(&ktreeMeanSigs[i * signatureSize], &means[node * signatureSize], signatureSize * sizeof(uint64_t));
			i++;
		}

		// get parents of ktree clusters
		set<size_t> parents;
		getNodesByLevel(root, 1, parents);
		size_t k = parents.size();
		vector<uint64_t> ktreeParentSigs(k * signatureSize);

		size_t j = 0;
		for (size_t node : parents) {
			memcpy(&ktreeParentSigs[j * signatureSize], &means[node * signatureSize], signatureSize * sizeof(uint64_t));
			j++;
		}


		// use the parent as seed, merge ktree 
		vector<uint64_t> meanSigs = ktreeParentSigs;
		vector<size_t> clusters(clusterCount);
		vector<vector<size_t>> clusterLists;
		for (size_t iteration = 0; iteration < 4; iteration++) {

			reclusterSignatures(clusters, meanSigs, ktreeMeanSigs, k);
			clusterLists = createClusterLists(clusters, k);
			meanSigs = createClusterSigs(clusterLists, ktreeMeanSigs, k);
		}

		auto it = parents.begin();
		for (size_t i = 0; i < clusterLists.size(); i++) {
			fprintf(stderr, "%zu,", *it);
			for (size_t node : clusterLists[i]) {
				fprintf(stderr, " %zu", node);
			}
			fprintf(stderr, "\n");
			it++;
			//RMSDs[i] = calcRMSD(&meanSigs[i*signatureSize], sigs, clusterLists[i]);
		}
		
		return clusters;
	}

};


void compressClusterList(vector<size_t> &clusters)
{
	unordered_map<size_t, size_t> remap;
	for (size_t &clus : clusters) {
		if (remap.count(clus)) {
			clus = remap[clus];
		}
		else {
			size_t newClus = remap.size();
			remap[clus] = newClus;
			clus = newClus;
		}
	}
	fprintf(stderr, "Output %zu clusters\n", remap.size());
}

vector<size_t> clusterSignatures(const vector<uint64_t> &sigs)
{
	size_t sigCount = sigs.size() / signatureSize;
	vector<size_t> clusters(sigCount);
	KTree tree(ktree_order, ktree_capacity);

	size_t firstNodes = 1;
	if (firstNodes > sigCount) firstNodes = sigCount;

	vector<size_t> insertionList;
	for (size_t i = 0; i <= firstNodes; i++) {
		insertionList.push_back(firstNodes - i);
	}

	default_random_engine rng;
	// Insert first 1 nodes single-threaded
	for (size_t i = 0; i < firstNodes; i++) {
		tree.insert(rng, &sigs[i * signatureSize], insertionList);
	}

	// What's the next free insertion point?
	size_t nextFree = insertionList.back();

	//#pragma omp parallel
	{
		default_random_engine rng;
		vector<size_t> insertionList;

		//#pragma omp for
		for (size_t i = nextFree; i < ktree_capacity; i++) {
			insertionList.push_back(ktree_capacity - i);
		}

		//#pragma omp for
		for (size_t i = firstNodes; i < sigCount; i++) {
			tree.insert(rng, &sigs[i * signatureSize], insertionList);
		}
	}

	// We've created the tree. Now reinsert everything
#pragma omp parallel for
	for (size_t i = 0; i < sigCount; i++) {
		size_t clus = tree.traverse(&sigs[i * signatureSize]);
		clusters[i] = clus;
	}

	vector<size_t>RMSDs = tree.updateTree(clusters, sigs);

	//// clustering the node centroids
	//vector<size_t> clustersOfClusters = tree.clusterClusters(rng, clusters, sigs);
	tree.restructureTree(rng, clusters, sigs);

	// We want to compress the cluster list down
	//compressClusterList(clusters);
	vector<size_t> compressedRMSDs = compressClusterRMSD(clusters, RMSDs);

	// output distortion
	string file_name = "silva-o" + to_string(ktree_order) +
		"-i0-distortion.txt";
	FILE * pFile = fopen(file_name.c_str(), "w");
	for (size_t i = 0; i < compressedRMSDs.size(); i++) {
		fprintf(pFile, "%zu,%zu\n", i, compressedRMSDs[i]);
	}

	// Recursively destroy all locks
	tree.destroyLocks();

	return clusters;

	//// We want to compress the cluster list down
	//compressClusterList(clusters);

	//// Recursively destroy all locks
	//tree.destroyLocks();

	//return clusters;
}

int main(int argc, char **argv)
{
	if (argc < 2) {
		fprintf(stderr, "Usage: %s (options) [fasta input]\n", argv[0]);
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -sw [signature width]\n");
		fprintf(stderr, "  -k [kmer length]\n");
		fprintf(stderr, "  -d [signature density]\n");
		fprintf(stderr, "  -o [tree order]\n");
		fprintf(stderr, "  -c [starting capacity]\n");
		fprintf(stderr, "  --fasta-output\n");
		return 1;
	}
	signatureWidth = 256;
	kmerLength = 5;
	density = 1.0f / 21.0f;
	fastaOutput = false;
	kmean_k = 100;

	string fastaFile = "";

	for (int a = 1; a < argc; a++) {
		string arg(argv[a]);
		if (arg == "-sw") signatureWidth = atoi(argv[++a]);
		else if (arg == "-k") kmerLength = atoi(argv[++a]);
		else if (arg == "-d") density = atof(argv[++a]);
		else if (arg == "-o") ktree_order = atoi(argv[++a]);
		else if (arg == "-c") ktree_capacity = atoi(argv[++a]);
		else if (arg == "-mk") kmean_k = atoi(argv[++a]);
		else if (arg == "--fasta-output") fastaOutput = true;
		else if (fastaFile.empty()) fastaFile = arg;
		else {
			fprintf(stderr, "Invalid or extra argument: %s\n", arg.c_str());
			exit(1);
		}
	}

	if (signatureWidth <= 0 || signatureWidth % 64 != 0) {
		fprintf(stderr, "Error: signature width is not a multiple of 64\n");
		return 1;
	}
	if (kmerLength <= 0) {
		fprintf(stderr, "Error: kmer length must be a positive nonzero integer\n");
		return 1;
	}
	if (kmean_k <= 0) {
		fprintf(stderr, "Error: kmean k-value must be a positive nonzero integer\n");
		return 1;
	}
	if (density < 0.0f || density > 1.0f) {
		fprintf(stderr, "Error: density must be a positive value between 0 and 1\n");
		return 1;
	}

	signatureSize = signatureWidth / 64;

	/*
	fprintf(stderr, "Loading fasta...");
	auto fasta = loadFasta(fastaFile.c_str());
	fprintf(stderr, " loaded %llu sequences\n", static_cast<unsigned long long>(fasta.size()));
	fprintf(stderr, "Converting fasta to signatures...");
	auto sigs = convertFastaToSignatures(fasta);
	fprintf(stderr, " done\n");
	*/

	fprintf(stderr, "Loading signatures...\n");
	auto sigs = readSignatures(fastaFile.c_str());

	//string file_name = "silva-o" + to_string(ktree_order) + "-i0.txt";
	//FILE * pFile = fopen(file_name.c_str(), "w");

	//fprintf(stderr, "Clustering signatures...\n");
	//auto clusters = clusterSignatures(sigs);
	//fprintf(stderr, "Writing output\n");
	//if (!fastaOutput) {
	//	outputClusters(pFile, clusters);
	//}
	///*else {
	//outputFastaClusters(clusters, fasta);
	//}*/

	//vector<size_t> orders = { 10,20,30,40,50,100,200,300 };
	vector<size_t> orders = { 10 };
	for (size_t order : orders) {
		ktree_order = order;
		string file_name = "silva-o" + to_string(ktree_order) + "-i0.txt";
		FILE * pFile = fopen(file_name.c_str(), "w");

		fprintf(stderr, "Clustering signatures...\n");
		auto clusters = clusterSignatures(sigs);
		fprintf(stderr, "Writing output\n");
		if (!fastaOutput) {
			outputClusters(pFile, clusters);
		}
	}




	return 0;
}