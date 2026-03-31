import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageDraw, ImageFont
import io

# --- 1. GLOBALE KONFIGURATION ---
st.set_page_config(page_title="Chemie-Designer Pro", page_icon="🧪", layout="wide")

# --- 2. FUNKTIONEN ---

def generate_mol_img(smiles):
    """Erzeugt ein Bild-Objekt mit fixierter Skalierung und dickeren Bindungen."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return None
    
    # R-Label Logik
    dummy_idx = 0
    labels = ["R", "R'", "R''", "R'''"]
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "*":
            label = labels[dummy_idx] if dummy_idx < len(labels) else f"R{dummy_idx}"
            atom.SetProp("atomLabel", label)
            dummy_idx += 1

    AllChem.Compute2DCoords(mol)
    
    # RDKit Drawing Options für "sauberen" Look
    d = rdMolDraw2D.MolDraw2DCairo(400, 400) # Etwas größer für Schärfe
    opts = d.drawOptions()
    opts.bondLineWidth = 3           # Dickere Bindungen
    opts.minFontSize = 20            # Lesbare Labels
    opts.fixedBondLength = 35        # Hält alle Moleküle im gleichen Maßstab
    
    d.DrawMolecule(mol)
    d.FinishDrawing()
    return Image.open(io.BytesIO(d.GetDrawingText()))

def export_to_rxn(smiles_str):
    """Erstellt einen MDL RXN Block für KingDraw/BioVIA."""
    try:
        rxn = AllChem.ReactionFromSmarts(smiles_str, useSmiles=True)
        return AllChem.ReactionToRxnBlock(rxn)
    except:
        return None

def parse_smiles(smiles_str):
    if not smiles_str or str(smiles_str).lower() == "nan": return None
    if ">>" in smiles_str:
        parts = smiles_str.split(">>")
        reactants = parts[0].split(".")
        products = parts[1].split(".")
        # Agent-Check (falls Format A.B>Agent>C)
        return {"r": reactants, "a": "", "p": products}
    return None

def draw_reaction_line(data, ratio_str, reaction_name="Unbenannte Reaktion"):
    r_parts = str(ratio_str).replace(">>", ".").split(".")
    canvas = Image.new('RGB', (2800, 600), color=(255, 255, 255))
    draw = ImageDraw.Draw(canvas)
    
    # Fonts laden mit massiv erhöhter Größe für Lesbarkeit
    try: 
        font_title = ImageFont.truetype("arial.ttf", 50)
        font_main = ImageFont.truetype("arial.ttf", 90) 
    except: 
        font_main = ImageFont.load_default()

    draw.text((50, 30), reaction_name, fill="black", font=font_title)
    
    x_cursor = 50
    y_midline = 300 # Das Zentrum der vertikalen Achse

    def draw_centered_text(text, x, y_mid, font, color="black"):
        bbox = draw.textbbox((0, 0), text, font=font)
        text_h = bbox[3] - bbox[1]
        draw.text((x, y_mid - (text_h // 2) - 10), text, fill=color, font=font)
        return bbox[2] - bbox[0] # Return width

    # 1. REAKTANDEN
    for i, sm in enumerate(data['r']):
        coeff = r_parts[i].strip() if i < len(r_parts) else ""
        if coeff not in ["1", "nan", "", "1.0"]:
            w = draw_centered_text(coeff, x_cursor, y_midline, font_main)
            x_cursor += w + 20
        
        mol_img = generate_mol_img(sm)
        if mol_img:
            mol_y = y_midline - (mol_img.height // 2)
            canvas.paste(mol_img, (x_cursor, mol_y), mol_img.convert("RGBA"))
            x_cursor += mol_img.width
        
        if i < len(data['r']) - 1:
            x_cursor += 20
            w = draw_centered_text("+", x_cursor, y_midline, font_main)
            x_cursor += w + 40

    # 2. PFEIL
    x_cursor += 40
    arrow_y = y_midline
    draw.line([(x_cursor, arrow_y), (x_cursor + 250, arrow_y)], fill="black", width=5)
    draw.polygon([(x_cursor + 250, arrow_y), (x_cursor + 230, arrow_y - 15), (x_cursor + 230, arrow_y + 15)], fill="black")
    x_cursor += 300

    # 3. PRODUKTE
    offset = len(data['r'])
    for i, sm in enumerate(data['p']):
        coeff_idx = offset + i
        coeff = r_parts[coeff_idx].strip() if coeff_idx < len(r_parts) else ""
        if coeff not in ["1", "nan", "", "1.0"]:
            w = draw_centered_text(coeff, x_cursor, y_midline, font_main)
            x_cursor += w + 20
        
        mol_img = generate_mol_img(sm)
        if mol_img:
            mol_y = y_midline - (mol_img.height // 2)
            canvas.paste(mol_img, (x_cursor, mol_y), mol_img.convert("RGBA"))
            x_cursor += mol_img.width
        
        if i < len(data['p']) - 1:
            x_cursor += 20
            w = draw_centered_text("+", x_cursor, y_midline, font_main)
            x_cursor += w + 40

    return canvas.crop((0, 0, x_cursor + 50, 600))

# --- 3. INTERFACE ---

st.title("⚗️ Chemie-Designer Pro")
tab1, tab2 = st.tabs(["✨ Einzel-Eingabe", "📂 Batch-Export"])

with tab1:
    col1, col2 = st.columns([2, 1])
    with col1:
        name = st.text_input("Reaktionsname", "Veresterung")
        smiles = st.text_input("SMILES (A.B>>C)", "C*(=O)O.OCC>>C*(=O)OCC")
        ratio = st.text_input("Verhältnisse (z.B. 2>>1)", "1>>1")

    if st.button("Generieren"):
        data = parse_smiles(smiles)
        if data:
            img = draw_reaction_line(data, ratio, name)
            st.image(img)
            
            # Downloads
            c1, c2 = st.columns(2)
            # Bild
            buf = io.BytesIO()
            img.save(buf, format="PNG")
            c1.download_button("💾 Bild speichern (PNG)", buf.getvalue(), f"{name}.png")
            
            # Chemie-Datei für KingDraw
            rxn_data = export_to_rxn(smiles)
            if rxn_data:
                c2.download_button("🧬 Editierbar (KingDraw/BioVIA)", rxn_data, f"{name}.rxn", "chemical/x-mdl-rxnfile")
        else:
            st.error("SMILES Format ungültig. Nutze '>>' als Trenner.")