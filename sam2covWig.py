#! /usr/bin/env python3

#__author__ = "Brandon Pickett"

# ----------- IMPORTS ---------------------------- ||
import sys
import argparse
import re
from collections import deque
from math import ceil

# ----------- CLASSES ---------------------------- ||
class SAMformatException(Exception):
	pass
	
# ---------- FUNCTIONS --------------------------- ||
def parseArgs():
	parser = argparse.ArgumentParser(prog="sam2covWig", description="Calculate coverage for SAM formatted-alignments and write them as a wiggle track with fixed steps", add_help=True, epilog="Incoming alignments should be provided via STDIN. Output wiggle file will be written to STDOUT. This program does not do any filtering; if filtering is desired, add appropriate options to your `samtools view` command. Alignments with a '*' in the CIGAR field will be ignored. All RNAMEs must be present in SAM @SQ headers with the LN tag defined. Currently, the M, D, X, and = CIGAR operations are used for coverage calculations. Note that SAM format (the input) uses a 1-based, fully-closed coordinate system: [start, end]. This is different from BED and PAF formats, which use a 0-based, half-open coordinate system: [start, end). WIG format (the output) uses the same coordinate system as SAM format.")
	parser.add_argument("-n", "--name", metavar="STR", type=str, action="store", dest="name", help="The name of this wiggle track, to be assigned as the value to the 'name' variable on the track definition line. If you supply the string 'Kitty cat', the track definition line will be 'track type=\"wiggle_0\" name=\"Kitty cat\"'. [default: Coverage]", default="Coverage", required=False)
	parser.add_argument("-w", "--window", metavar="INT", type=int, action="store", dest="window", help="The size of the window for calculating coverage. This value will be used while making the coverage calculations, and it will be supplied as the value to the 'step' and 'span' variables in the step definitions lines. If you supply a value of 128, the step definition lines will be 'fixedStep chrom=<chrC> step=128 span=128', where '<chrC>' is replaced by the RNAME (column 3 (1-based)) value for the current SAM alignment record, e.g., chr1, chr2, chrX, etc. [default: 1024]", default=1024, required=False)
	args = parser.parse_args()
	return args

def sumTargetConsumingCigarOps(cigar):
	TARGET_CONSUMERS_REGEX = r'([0-9]+)[=DMX]'
	consumed_target_bases = 0
	for match in re.finditer(TARGET_CONSUMERS_REGEX, cigar):
		consumed_target_bases += int(match.group(1))
	return consumed_target_bases

def incrementCountsForCoveredWindows(windows, current_window_start, target_start, target_end, target_len, window_size):
	# windows is modified directly (w/o needing to be returned).
	# the other variables are numbers and are neither modified nor returned.
	# The current window (which may have to be added to the windows list) will be
	#	incremented along with subsequent covered windows.

	#assert target_len > 0
	#assert current_window_start > 0
	#assert window_size > 0
	#assert target_start <= current_window_start + window_size
	#assert target_start >= current_window_start
	#assert target_end == target_start + target_len - 1

	max_possible_windows_covered = int(target_len / window_size) + 1
	num_windows_covered = ceil((target_end + 1 - current_window_start) / window_size)
	#num_windows_covered = sum(1 for i in range(current_window_start, target_end + 1, window_size))

	# insert extra windows into list if needed
	while len(windows) < num_windows_covered:
		windows.append(0)
	
	# "loop" through covered windows to increment the count
	windows.rotate(-num_windows_covered)
	for i in range(num_windows_covered):
		windows.appendleft(windows.pop() + 1)


# ------------- MAIN ----------------------------- ||
if __name__ == "__main__":
	
	# handle the arguments
	args = parseArgs()
	WINDOW_SIZE = args.window
	TRACK_NAME = args.name

	# write the wig track definition line
	print(f'track type="wiggle_0" name="{TRACK_NAME}"', file=sys.stdout)

	# process the input SAM
	line = sys.stdin.readline()

	# 	process header lines
	SQ_SN_REGEX = r"([LS]N):([^\t]+)"
	chr_lens = {}
	while line[0] == '@':
		if line.startswith("@SQ"):
			temp_pair = dict(re.findall(SQ_SN_REGEX, line))
			if not ("SN" in temp_pair and "LN" in temp_pair): raise SAMformatException("SAM @SQ header record was missing SN and/or LN")
			chr_lens[temp_pair["SN"]] = int(temp_pair["LN"])
		line = sys.stdin.readline()
	
	# 	process the tabular lines
	windows = deque() # for use as a linked-list, inserted values are counts of reads that cover at least one base in that window
	current_chr = ''
	current_chr_len = 0
	current_window_start = 1 # 1-based
	while line != '':
		fields = line.rstrip('\n').split('\t')

		# extract needed info
		target_name = fields[2]
		target_start = int(fields[3]) # 1-based start position
		cigar = fields[5]

		# skip cigar == '*' lines (presumably because they don't map; if all cigars are '*', the coverages will be 0 in all windows)
		if cigar != '*':

			# if we've moved on to another chromosome (i.e., target sequence)
			if current_chr != target_name:
				# finish previous chr's windows
				#	advance current windows start position based on num stored windows
				current_window_start += len(windows) * WINDOW_SIZE
				#	write each window's current value
				for window in windows:
					print(window, file=sys.stdout)
				#	empty the windows list
				windows.clear()
				#	write unstored windows' values (all 0)
				for i in range(current_window_start, current_chr_len + 1, WINDOW_SIZE):
					print(0, file=sys.stdout)

				# set things up for this new chr
				# 	set current chr
				current_chr = target_name
				#	set current chr len
				if not current_chr in chr_lens: raise SAMformatException("SAM alignment record had RNAME value not present in SAM @SQ headers")
				current_chr_len = chr_lens[current_chr]
				# 	reset window start
				current_window_start = 1 # 1-based, fully-closed: [start, end]
				# 	write wig step declaration line
				print(f"fixedStep chrom={current_chr} start={current_window_start} step={WINDOW_SIZE} span={WINDOW_SIZE}", file=sys.stdout)

			# advance current window to this SAM alignment's start position, if we aren't already in the correct window
			while target_start > current_window_start + WINDOW_SIZE - 1:
				# current SAM alignment starts outside (to the right of) the current window. So...
				current_window_coverage = 0

				# if we have any stored coverages for windows
				if len(windows): # if len(windows) > 0
					# extract current window's coverage while removing window from list of windows
					current_window_coverage = windows.popleft()

				# write the window (even if it is zero). Note that we do this regardless of whether
				#	we had a value stored for the window in the windows list (which means the current
				#	window's coverage is 0). This is because we may have a coverage gap where no
				#	counts were recorded.
				print(current_window_coverage, file=sys.stdout)

				# advance the current_window_start by the step size. Note that we do this regardless of
				#	whether we have stored values in the windows list. This is because we may have a
				#	coverage gap where no counts were recorded.
				current_window_start += WINDOW_SIZE

			# calculate length of alignment on target from cigar
			target_len = sumTargetConsumingCigarOps(cigar)

			# calculate target end position
			target_end = target_start + target_len -1 # 1-based, fully-closed: [start, end]
			#target_end = target_start + target_len # this would be for half open style: [start, end)

			# increment counts for current window (covered by definition) and subsequent covered windows
			incrementCountsForCoveredWindows(windows, current_window_start, target_start, target_end, target_len, WINDOW_SIZE)

		# advance to next line
		line = sys.stdin.readline()
	
	# write remaining windows (if any) now that we have processed all SAM records
	#	advance current windows start position based on num stored windows
	current_window_start += len(windows) * WINDOW_SIZE
	#	write each window's current value
	for window in windows:
		print(window, file=sys.stdout)
	#	empty the windows list
	windows.clear()
	#	write unstored windows' values (all 0)
	for i in range(current_window_start, current_chr_len + 1, WINDOW_SIZE):
		print(0, file=sys.stdout)

