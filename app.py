import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageDraw, ImageFont
import io
import os

# --- 1. GLOBALE KONFIGURATION (MUSS ZUERST KOMMEN) ---
st.set_page_config(page_title="Chemie-Designer", page_icon="🧪", layout="wide")

# --- 2. FUNKTIONEN ---

def generate_mol_img(smiles):
    """Erzeugt ein Bild-Objekt aus einem SMILES und ersetzt '*' durch R-Labels."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    
    dummy_idx = 0
    labels = ["R", "R'", "R''", "R'''"]
    
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "*":
            label = labels[dummy_idx] if dummy_idx < len(labels) else f"R{dummy_idx}"
            atom.SetProp("atomLabel", label)
            dummy_idx += 1

    AllChem.Compute2DCoords(mol)
    
    d = rdMolDraw2D.MolDraw2DCairo(300, 300)
    opts = d.drawOptions()
    opts.bondLineWidth = 2
    opts.prepareMolsBeforeDrawing = False 
    d.DrawMolecule(mol)
    d.FinishDrawing()
    
    return Image.open(io.BytesIO(d.GetDrawingText()))

def parse_smiles(smiles_str):
    if not smiles_str or str(smiles_str).lower() == "nan" or str(smiles_str).strip() == "":
        return None
    parts = str(smiles_str).split(">")
    if len(parts) == 3:
        return {"r": parts[0].split("."), "a": parts[1], "p": parts[2].split(".")}
    elif ">>" in smiles_str:
        temp = smiles_str.split(">>")
        return {"r": temp[0].split("."), "a": "", "p": temp[1].split(".")}
    return None

def draw_reaction_line(data, ratio_str, reaction_name="Unbenannte Reaktion"):
    """Erstellt das kombinierte Bild mit Titel oben drauf."""
    r_parts = str(ratio_str).replace(">>", ".").split(".")
    
    canvas = Image.new('RGB', (2500, 500), color=(255, 255, 255))
    draw = ImageDraw.Draw(canvas)
    
    try: 
        font_title = ImageFont.truetype("arial.ttf", 45)
        font_main = ImageFont.truetype("arial.ttf", 35)
    except: 
        font_title = ImageFont.load_default()
        font_main = ImageFont.load_default()

    # 0. REAKTIONSNAME ZEICHNEN
    draw.text((40, 20), reaction_name, fill="black", font=font_title)
    
    x_cursor, y_center = 40, 150 

    # 1. REAKTANDEN
    for i, sm in enumerate(data['r']):
        coeff = r_parts[i].strip() if i < len(r_parts) else ""
        if coeff not in ["1", "nan", "", "1.0"]:
            draw.text((x_cursor, y_center + 130), coeff, fill="black", font=font_main)
            x_cursor += 50
        
        mol_img = generate_mol_img(sm)
        if mol_img:
            canvas.paste(mol_img, (x_cursor, y_center), mol_img.convert("RGBA"))
            x_cursor += 310
        
        if i < len(data['r']) - 1:
            draw.text((x_cursor + 10, y_center + 130), "+", fill="black", font=font_main)
            x_cursor += 80

    # 2. PFEIL
    x_cursor += 20
    arrow_len = 220
    y_f = y_center + 150
    draw.line([(x_cursor, y_f), (x_cursor + arrow_len, y_f)], fill="black", width=4)
    draw.polygon([(x_cursor + arrow_len, y_f), (x_cursor + arrow_len - 15, y_f - 8), (x_cursor + arrow_len - 15, y_f + 8)], fill="black")
    if data['a']:
        draw.text((x_cursor + 10, y_f - 60), str(data['a']), fill="blue", font=font_main)
    x_cursor += arrow_len + 40

    # 3. PRODUKTE
    offset = len(data['r'])
    for i, sm in enumerate(data['p']):
        coeff = r_parts[offset + i].strip() if (offset + i) < len(r_parts) else ""
        if coeff not in ["1", "nan", "", "1.0"]:
            draw.text((x_cursor, y_center + 130), coeff, fill="black", font=font_main)
            x_cursor += 50
        
        mol_img = generate_mol_img(sm)
        if mol_img:
            canvas.paste(mol_img, (x_cursor, y_center), mol_img.convert("RGBA"))
            x_cursor += 310
        
        if i < len(data['p']) - 1:
            draw.text((x_cursor + 10, y_center + 130), "+", fill="black", font=font_main)
            x_cursor += 80

    return canvas.crop((0, 0, x_cursor + 20, 500))

# --- 3. HAUPT-INTERFACE ---

st.markdown("# ⚗️ Chemie-Reaktions Visualisierer")
st.write("Erstelle saubere Reaktionsgleichungen mit Namen und R-Gruppen.")
st.divider()

# Tabs definieren
tab1, tab2 = st.tabs(["✨ Einzelne Eingabe", "📂 Excel-Verarbeitung"])

# --- TAB 1: EINZEL-EINGABE ---
with tab1:
    col1, col2 = st.columns([2, 1])
    with col1:
        reaction_name_input = st.text_input("Name der Reaktion", "Meine Reaktion")
        smiles_input = st.text_input("SMILES String (z.B. A.B>>C)", "C*(=O)O.OCC>>C*(=O)OCC")
        ratio_input = st.text_input("Molare Verhältnisse", "1>>1")
    
    if st.button("Visualisieren", key="btn_single"):
        data = parse_smiles(smiles_input)
        if data:
            result_img = draw_reaction_line(data, ratio_input, reaction_name_input)
            st.image(result_img, caption="Vorschau")
            
            # Download vorbereiten
            buf = io.BytesIO()
            result_img.save(buf, format="PNG")
            filename = f"{reaction_name_input.replace(' ', '_')}.png"
            st.download_button("Bild herunterladen", buf.getvalue(), filename, "image/png")
        else:
            st.error("Ungültiges SMILES-Format. Bitte prüfe deine Eingabe.")

# --- TAB 2: EXCEL-VERARBEITUNG ---
with tab2:
    uploaded_file = st.file_uploader("Excel-Datei hochladen (.xlsx)", type=["xlsx"])
    if uploaded_file:
        df = pd.read_excel(uploaded_file)
        st.write("Vorschau der geladenen Daten:")
        st.dataframe(df.head())
        
        if st.button("Alle Reaktionen generieren", key="btn_batch"):
            for idx, row in df.iterrows():
                nom = str(row.get('Nom', f"Reaction_{idx}"))
                ratio = str(row.get('Ratio', '1>>1'))
                smiles_val = row.get('SMILES_Allg')
                
                data = parse_smiles(smiles_val)
                if data:
                    st.write(f"### {nom}")
                    img = draw_reaction_line(data, ratio, nom)
                    st.image(img)
                    
                    # Optional: Kleiner Download-Link pro Bild
                    buf = io.BytesIO()
                    img.save(buf, format="PNG")
                    st.download_button(f"Download {nom}", buf.getvalue(), f"{nom}.png", "image/png", key=f"dl_{idx}")
                else:
                    st.warning(f"Konnte SMILES in Zeile {idx} ({nom}) nicht lesen.")