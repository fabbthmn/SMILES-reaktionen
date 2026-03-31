import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageDraw, ImageFont
import io

# --- 1. GLOBALE KONFIGURATION ---
st.set_page_config(page_title="Chemie-Designer Pro", page_icon="🧪", layout="wide")

# --- 2. KERN-FUNKTIONEN ---

def generate_mol_img(smiles):
    """Generiert ein einzelnes Molekülbild mit RDKit."""
    if not smiles or str(smiles).lower() == "nan": return None
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
    
    d = rdMolDraw2D.MolDraw2DCairo(550, 550) 
    opts = d.drawOptions()
    opts.bondLineWidth = 5
    opts.minFontSize = 35
    opts.fixedBondLength = 50       
    
    d.DrawMolecule(mol)
    d.FinishDrawing()
    return Image.open(io.BytesIO(d.GetDrawingText()))

def parse_smiles(smiles_str):
    """Trennt SMILES in Reaktanten und Produkte."""
    if not smiles_str or str(smiles_str).lower() == "nan": return None
    if ">>" in smiles_str:
        parts = smiles_str.split(">>")
        return {"r": parts[0].split("."), "p": parts[1].split(".")}
    return None

def draw_single_reaction_line(data, ratio_str, font_main, font_title=None, reaction_name=None):
    """Zeichnet eine einzelne Zeile (eine SMILES-Variante)."""
    r_parts = str(ratio_str).replace(">>", ".").split(".")
    
    anzahl_elemente = len(data['r']) + len(data['p'])
    dynamische_breite = (anzahl_elemente * 650) + 1200 
    dynamische_breite = max(dynamische_breite, 1800) 
    
    line_h = 850 
    canvas = Image.new('RGB', (dynamische_breite, line_h), color=(255, 255, 255))
    draw = ImageDraw.Draw(canvas)
    
    x_cursor = 80
    y_midline = line_h // 2

    if reaction_name and font_title:
        draw.text((80, 40), str(reaction_name), fill="black", font=font_title)

    def draw_centered_text(text, x, y_mid, font):
        bbox = draw.textbbox((0, 0), text, font=font)
        text_w = bbox[2] - bbox[0]
        text_h = bbox[3] - bbox[1]
        draw.text((x, y_mid - (text_h // 2) - 25), text, fill="black", font=font)
        return text_w

    # 1. REAKTANDEN
    for i, sm in enumerate(data['r']):
        coeff = r_parts[i].strip() if i < len(r_parts) else ""
        if coeff not in ["1", "nan", "", "1.0"]:
            w = draw_centered_text(coeff, x_cursor, y_midline, font_main)
            x_cursor += w + 50 
        
        mol_img = generate_mol_img(sm)
        if mol_img:
            mol_y = y_midline - (mol_img.height // 2)
            canvas.paste(mol_img, (x_cursor, mol_y), mol_img.convert("RGBA"))
            x_cursor += mol_img.width
        
        if i < len(data['r']) - 1:
            x_cursor += 50
            w = draw_centered_text("+", x_cursor, y_midline, font_main)
            x_cursor += w + 70

    # 2. PFEIL
    x_cursor += 70
    draw.line([(x_cursor, y_midline), (x_cursor + 350, y_midline)], fill="black", width=10)
    draw.polygon([(x_cursor + 350, y_midline), (x_cursor + 300, y_midline - 30), (x_cursor + 300, y_midline + 30)], fill="black")
    x_cursor += 450

    # 3. PRODUKTE
    offset = len(data['r'])
    for i, sm in enumerate(data['p']):
        coeff_idx = offset + i
        coeff = r_parts[coeff_idx].strip() if coeff_idx < len(r_parts) else ""
        if coeff not in ["1", "nan", "", "1.0"]:
            w = draw_centered_text(coeff, x_cursor, y_midline, font_main)
            x_cursor += w + 50
        
        mol_img = generate_mol_img(sm)
        if mol_img:
            mol_y = y_midline - (mol_img.height // 2)
            canvas.paste(mol_img, (x_cursor, mol_y), mol_img.convert("RGBA"))
            x_cursor += mol_img.width
        
        if i < len(data['p']) - 1:
            x_cursor += 50
            w = draw_centered_text("+", x_cursor, y_midline, font_main)
            x_cursor += w + 70

    return canvas.crop((0, 0, x_cursor + 120, line_h))

def combine_reactions(smiles_list, ratio_list, reaction_name):
    """Kombiniert mehrere Varianten vertikal in ein Gesamtbild."""
    try: 
        font_title = ImageFont.truetype("arial.ttf", 80)
        font_main = ImageFont.truetype("arial.ttf", 100)
    except: 
        font_title = ImageFont.load_default()
        font_main = ImageFont.load_default()

    images = []
    for i, sm in enumerate(smiles_list):
        data = parse_smiles(sm)
        if data:
            n = reaction_name if i == 0 else None
            images.append(draw_single_reaction_line(data, ratio_list[i], font_main, font_title, n))
    
    if not images: return None

    max_w = max(img.width for img in images)
    total_h = sum(img.height for img in images) + 60
    
    combined = Image.new('RGB', (max_w, total_h), (255, 255, 255))
    y_offset = 0
    draw = ImageDraw.Draw(combined)
    
    for i, img in enumerate(images):
        combined.paste(img, (0, y_offset))
        y_offset += img.height
        if i < len(images) - 1:
            # HIER WAR DER FEHLER: Korrekte Klammersetzung für draw.line
            draw.line([(80, y_offset), (max_w - 80, y_offset)], fill=(210, 210, 210), width=4)
            y_offset += 30

    return combined

# --- 3. STREAMLIT INTERFACE ---

st.title("⚗️ Chemie-Designer Pro (Multi-SMILES)")
tab1, tab2 = st.tabs(["✨ Einzel-Eingabe", "📂 Batch-Export"])

with tab1:
    with st.expander("📝 Eingabe-Details", expanded=True):
        name = st.text_input("Reaktionsname", "Beispielreaktion")
        st.markdown("### SMILES-Varianten")
        c1, c2 = st.columns(2)
        s1 = c1.text_input("SMILES Variante 1", "C*(=O)O.OCC>>C*(=O)OCC")
        r1 = c1.text_input("Verhältnis 1", "1>>1")
        s2 = c2.text_input("SMILES Variante 2 (Optional)", "")
        r2 = c2.text_input("Verhältnis 2", "1>>1")

    if st.button("🚀 Generieren", use_container_width=True):
        s_list = [s for s in [s1, s2] if s.strip()]
        r_list = [r1, r2]
        img = combine_reactions(s_list, r_list, name)
        if img:
            st.image(img)
            buf = io.BytesIO()
            img.save(buf, format="PNG")
            st.download_button("💾 Bild herunterladen", buf.getvalue(), f"{name}.png")

with tab2:
    uploaded_file = st.file_uploader("Excel-Datei (.xlsx) hochladen", type=["xlsx"])
    if uploaded_file:
        df = pd.read_excel(uploaded_file)
        cols = df.columns.tolist()
        c_a, c_b, c_c, c_d = st.columns(4)
        col_name = c_a.selectbox("Spalte: Name", cols)
        col_s1 = c_b.selectbox("Spalte: SMILES A", cols)
        col_s2 = c_c.selectbox("Spalte: SMILES B (Optional)", ["Keine"] + cols)
        col_ratio = c_d.selectbox("Spalte: Verhältnisse", ["Keine"] + cols)

        if st.button("🏁 Alle Reaktionen verarbeiten"):
            for idx, row in df.iterrows():
                nom = str(row[col_name])
                s_list = [str(row[col_s1])]
                r_val = str(row[col_ratio]) if col_ratio != "Keine" else "1>>1"
                if col_s2 != "Keine" and pd.notna(row[col_s2]):
                    s_list.append(str(row[col_s2]))
                
                img = combine_reactions(s_list, [r_val, r_val], nom)
                if img:
                    st.markdown("---")
                    st.subheader(nom)
                    st.image(img)
                    buf = io.BytesIO()
                    img.save(buf, format="PNG")
                    st.download_button(f"Download {nom}", buf.getvalue(), f"{nom}.png", key=f"dl_{idx}")