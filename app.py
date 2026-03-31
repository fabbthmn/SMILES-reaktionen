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
    if not smiles or str(smiles).lower() == "nan": return None
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return None
    
    dummy_idx = 0
    labels = ["R", "R'", "R''", "R'''"]
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "*":
            label = labels[dummy_idx] if dummy_idx < len(labels) else f"R{dummy_idx}"
            atom.SetProp("atomLabel", label)
            dummy_idx += 1

    AllChem.Compute2DCoords(mol)
    d = rdMolDraw2D.MolDraw2DCairo(500, 500)
    opts = d.drawOptions()
    opts.bondLineWidth = 4
    opts.minFontSize = 25
    opts.fixedBondLength = 45
    
    d.DrawMolecule(mol)
    d.FinishDrawing()
    return Image.open(io.BytesIO(d.GetDrawingText()))

def parse_smiles(smiles_str):
    if not smiles_str or str(smiles_str).lower() == "nan": return None
    if ">>" in smiles_str:
        parts = smiles_str.split(">>")
        reactants = parts[0].split(".")
        products = parts[1].split(".")
        return {"r": reactants, "p": products}
    return None

def draw_single_reaction_line(data, ratio_str, font_main, font_title=None, reaction_name=None):
    """Hilfsfunktion: Zeichnet EINE Zeile einer Reaktion."""
    r_parts = str(ratio_str).replace(">>", ".").split(".")
    
    # Leinwand-Setup
    anzahl_elemente = len(data['r']) + len(data['p'])
    dynamische_breite = (anzahl_elemente * 550) + 800 
    dynamische_breite = max(dynamische_breite, 1500) 
    
    line_h = 600 # Höhe einer einzelnen Zeile
    canvas = Image.new('RGB', (dynamische_breite, line_h), color=(255, 255, 255))
    draw = ImageDraw.Draw(canvas)
    
    x_cursor = 50
    y_midline = line_h // 2

    if reaction_name and font_title:
        draw.text((50, 20), str(reaction_name), fill="black", font=font_title)

    def draw_centered_text(text, x, y_mid, font):
        bbox = draw.textbbox((0, 0), text, font=font)
        text_h = bbox[3] - bbox[1]
        draw.text((x, y_mid - (text_h // 2)), text, fill="black", font=font)
        return bbox[2] - bbox[0]

    # Reaktanten
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

    # Pfeil
    x_cursor += 40
    draw.line([(x_cursor, y_midline), (x_cursor + 250, y_midline)], fill="black", width=5)
    draw.polygon([(x_cursor + 250, y_midline), (x_cursor + 230, y_midline - 15), (x_cursor + 230, y_midline + 15)], fill="black")
    x_cursor += 300

    # Produkte
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

    return canvas.crop((0, 0, x_cursor + 50, line_h))

def combine_reactions(smiles_list, ratio_list, reaction_name):
    """Kombiniert mehrere SMILES-Varianten in ein Bild untereinander."""
    # Schriften laden
    try: 
        font_title = ImageFont.truetype("arial.ttf", 45)
        font_main = ImageFont.truetype("arial.ttf", 65)
    except: 
        font_title = ImageFont.load_default()
        font_main = ImageFont.load_default()

    images = []
    for i, sm in enumerate(smiles_list):
        data = parse_smiles(sm)
        if data:
            # Nur beim ersten Bild den Namen oben einfügen
            name_to_draw = reaction_name if i == 0 else None
            images.append(draw_single_reaction_line(data, ratio_list[i], font_main, font_title, name_to_draw))
    
    if not images: return None

    # Kombiniere Bilder vertikal
    max_w = max(img.width for img in images)
    total_h = sum(img.height for img in images) + 20
    
    combined = Image.new('RGB', (max_w, total_h), (255, 255, 255))
    y_offset = 0
    for img in images:
        combined.paste(img, (0, y_offset))
        y_offset += img.height
        # Kleine Trennlinie zwischen Varianten
        draw = ImageDraw.Draw(combined)
        draw.line([(0, y_offset), (max_w, y_offset)], fill=(220, 220, 220), width=2)
        y_offset += 5

    return combined

# --- 3. INTERFACE ---
st.title("⚗️ Chemie-Designer Pro (Multi-SMILES)")
tab1, tab2 = st.tabs(["✨ Einzel-Eingabe", "📂 Batch-Export"])

with tab1:
    with st.expander("Eingaben", expanded=True):
        col1, col2 = st.columns(2)
        name = col1.text_input("Reaktionsname", "Veresterung")
        
        st.markdown("---")
        c1, c2 = st.columns(2)
        smiles1 = c1.text_input("SMILES Variante A", "C*(=O)O.OCC>>C*(=O)OCC")
        ratio1 = c1.text_input("Verhältnis A", "1>>1")
        
        smiles2 = c2.text_input("SMILES Variante B (Optional)", "")
        ratio2 = c2.text_input("Verhältnis B", "1>>1")

    if st.button("Generieren"):
        s_list = [s for s in [smiles1, smiles2] if s]
        r_list = [ratio1, ratio2]
        
        img = combine_reactions(s_list, r_list, name)
        if img:
            st.image(img)
            buf = io.BytesIO()
            img.save(buf, format="PNG")
            st.download_button("💾 Kombiniertes Bild speichern", buf.getvalue(), f"{name}_multi.png")

with tab2:
    uploaded_file = st.file_uploader("Excel-Datei hochladen", type=["xlsx"])
    if uploaded_file:
        df = pd.read_excel(uploaded_file)
        cols = df.columns.tolist()
        
        c1, c2, c3, c4 = st.columns(4)
        col_name = c1.selectbox("Spalte Name", cols)
        col_s1 = c2.selectbox("Spalte SMILES A", cols)
        col_s2 = c3.selectbox("Spalte SMILES B (Optional)", ["Keine"] + cols)
        col_ratio = c4.selectbox("Spalte Verhältnis", ["Keine"] + cols)

        if st.button("🚀 Batch Prozess starten"):
            for idx, row in df.iterrows():
                nom = str(row[col_name])
                s_list = [str(row[col_s1])]
                r_list = [str(row[col_ratio]) if col_ratio != "Keine" else "1>>1"]
                
                if col_s2 != "Keine" and pd.notna(row[col_s2]):
                    s_list.append(str(row[col_s2]))
                    r_list.append(r_list[0]) # Nutzt gleiches Verhältnis oder Logik erweitern
                
                img = combine_reactions(s_list, r_list, nom)
                if img:
                    st.subheader(f"Reaktion: {nom}")
                    st.image(img)
                    # Download...