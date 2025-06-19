
# Loading FASTA file and extracting target amplicon based user preferences

from Bio import SeqIO

def handle_fasta(fasta_path):
    
    # Loading DNA Sequence from FASTA file
   
    try:
        record = next(SeqIO.parse(fasta_path, "fasta"))
        #print(record)
        gene_seq = str(record.seq).upper()
        #print(gene_seq)
        print(f"\nLoaded sequence of length {len(gene_seq)} from {record.id}")
   
    except Exception as e:
        print(f"\nError reading FASTA file: {e}")
        return None, None, None
    
    return gene_seq


def extract_target(gene_seq, min_target_size=100, max_target_size=200, flank=50):
    # Ask user for mode of target region selection
    # Validating target selection criteria
    # Generating the target sequence

    while True: 

        print("\nSelect target extraction method:")
        print("1. By position (start-end)")
        print("2. By sequence match")
        choice = int(input("Enter 1 or 2: ").strip())
        
        if choice == 1:
            while True:
                start = int(input("\nEnter target start position: ").strip()) - 1
                end = int(input("Enter target end position: ").strip()) - 1
                target_size = end - start + 1

                if start < 0 or end > (len(gene_seq)-1): # Validating co-ordinate position
                    print("\nEntered coordinates exceed the bounds of the sequence. Please enter valid coordinates.")
                    continue
                if start >= end: # Validating if start position comes before end position
                    print("\nStart position cannot be greater than the end position. Please enter valid coordinates.")
                    continue
                if target_size not in range(min_target_size, max_target_size+1): # Validating target size
                    print(f"\nYou have selected a target region of {target_size} bp.")
                    print(f"Target size must be in the range from {min_target_size} to {max_target_size} bp. Please enter valid coordinates.")
                    continue
                if (start-flank-1) < 0 or (end+flank) > (len(gene_seq)-1): # Validating presence of sufficent flanking sequence
                    print(f"Insufficient flanking sequence: at least {flank} bp is required on both ends of the selected target region for primer binding. Please enter valid coordinates.")
                    continue
                else:
                    print(f"\nYou have selected a target region of {target_size} bp.")
                    target_seq = gene_seq[start:end+1]
                    print(f"Selected targeted sequence:\n{target_seq}")
                    return target_seq, start, end, flank
        
        elif choice == 2:
            while True:
                target_seq = str(input("\nEnter target sequence: ")).strip().upper()
                target_size = len(target_seq)
                idx = gene_seq.find(target_seq)
                
                if idx == -1:
                    print("Target region not found in sequence.")
                    continue
                else:
                    start = idx
                    end = idx + len(target_seq) - 1
                    if (end-start) not in range(min_target_size, max_target_size+1):
                        print(f"\nYou have selected a target region of {target_size} bp.")
                        print(f"Target size must be in the range from {min_target_size} to {max_target_size} bp. Please enter valid target sequence.")
                        continue
                    if (start-flank-1) < 0 or (end+flank) > (len(gene_seq)-1):
                        print(f"Insufficient flanking sequence: at least {flank} bp is required on both ends of the selected target region for primer binding. Please enter valid target sequence.")
                        continue
                    else:
                        print(f"\nYou have selected a target region of {target_size} bp, at positions {start+1}-{end+1}.")
                        return target_seq, start, end, flank

        else:
            print("Invalid choice.")
