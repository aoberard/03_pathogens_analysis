import sys
import pandas as pd
from pathlib import Path
import os
import argparse

    # if len(sys.argv) < 3:
    #     print("Usage : python nom_du_script.py <chemin_fichier_input (obligatoire)> <separateur fichier = ',' ou ';' ou '\\t' (obligatoire)> <taxon_target (obligatoire)> <prefixe (facultatif)>")
    #     sys.exit(1)
        
    # elif len(sys.argv)>4:
    #     prefixe = sys.argv[4]
        
    # print(len(sys.argv))
    
    # Récupération des arguments
#     input_file = sys.argv[1]
#     separateur = sys.argv[2]
#     taxon_target = sys.argv[3]

def main():
    # Vérification des arguments
    parser=argparse.ArgumentParser()
    parser.add_argument("--input", help="fichier input csv, txt ou excel (OBLIGATOIRE)")
    parser.add_argument("--sep", help="separateur utilisée dans le fichier parmis : ',' ou ';' ou '\\t' (OBLIGATOIRE)")
    parser.add_argument("--target", help="taxon cible (OBLIGATOIRE)")
    parser.add_argument("--pref", help="prefixe a utiliser dans l'en-tete fasta (FACULTATIF)")

    # print(parser.format_help())
    args=parser.parse_args()



    input_file = args.input
    separateur = args.sep
    taxon_target = args.target
    prefixe = args.pref

    # Vérification de l'existence du fichier Excel
    excel_path = Path(input_file)
    if not excel_path.exists():
        print(f"Erreur : Le fichier {input_file} n'existe pas.")
        sys.exit(1)

    # Lecture du fichier Excel avec gestion des engines
    try:
        if (input_file.endswith(".txt") or input_file.endswith(".csv")) and (separateur in [',' , ';', '\t']):
            df = pd.read_csv(excel_path, sep=f"{separateur}") 
        elif input_file.endswith(".xlsx"):
            df = pd.read_excel(excel_path, engine="openpyxl")  # Pour les fichiers modernes
        elif input_file.endswith(".xls"):
            df = pd.read_excel(excel_path, engine="xlrd")  # Pour les anciens fichiers Excel
        else:
            print("Erreur : Format de fichier non supporté ou pas de séparateur renseigné.")
            sys.exit(1)
    except Exception as e:
        print(f"Erreur lors de la lecture du fichier Excel : {e}")
        sys.exit(1)

    # Vérification des colonnes nécessaires
    required_columns = {"genus", "seed_sequence", "observation_name"}
    if not required_columns.issubset(df.columns):
        print(f"Erreur : Le fichier Excel doit contenir les colonnes suivantes : {', '.join(required_columns)}")
        sys.exit(1)

    # Filtrage des données pour le taxon cible
    filtered_df = df[(df["genus"]==taxon_target) & ((df["species"].str.contains("Multi-affiliation")) | (df["species"].str.contains("unknown"))) & (df["seed_sequence"].notna())]
    
    if filtered_df.empty:
        print(f"Erreur : aucunes séquences trouvées où l'espèce est unknown ou Multi-affiliation {taxon_target}.")
        sys.exit(1)

    # Création du fichier FASTA
    output_fasta = Path.cwd() / f"{taxon_target}_sequences.fasta"
    
    if os.path.isfile(output_fasta):
        print(f"Erreur : traitement annulé, attention le fichier fasta {output_fasta} existe déjà.")
        sys.exit(1)
        
    with open(output_fasta, "w") as fasta_file:
        for _, row in filtered_df.iterrows():
            if prefixe:
                header = f">{prefixe}_{row['genus']}_{row['observation_name']}"
            else:
                header = f">{row['genus']}_{row['observation_name']}"
            sequence = row["seed_sequence"]
            fasta_file.write(f"{header}\n{sequence}\n")

    # Résumé
    total_rows = len(df)
    generated_sequences = len(filtered_df)
    nb_unknown = len(filtered_df[filtered_df["species"].str.contains("unknown")])
    nb_multi = len(filtered_df[filtered_df["species"].str.contains("Multi-affiliation")])

    print(f"Résumé :")
    print(f"- Nombre de lignes dans le tableau de départ : {total_rows}")
    print(f"- Nombre de séquences générées dans le fichier fasta : {generated_sequences}")
    print(f"- Dont nombre de séquences unknown species : {nb_unknown}")
    print(f"- Dont nombre de séquences Multi-affiliation : {nb_multi}")
    print(f"Fichier FASTA créé : {output_fasta}")

if __name__ == "__main__":
    main()