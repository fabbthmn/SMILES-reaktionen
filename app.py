import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageDraw, ImageFont
import io
import zipfile
import os
import tempfile

# --- 1. GLOBALE KONFIGURATION ---
st.set_page_config(page_title="Chemie-Designer Pro", page_icon="🧪", layout="wide")

# --- 2. HILFSFUNKTIONEN ---

@st.cache_resource
def get_font(size):
    """Lädt die Schriftart und speichert sie im Cache (@st.cache_resource für UI-Elemente)."""
    font_paths = [
        "arial.ttf", 
        "OpenSans-Regular.ttf", # Besser für Cloud-Deployments (Datei muss im Ordner liegen)
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 
        "/System/Library/Fonts/Helvetica.ttc"
    ]
    for path in font_paths:
        try:
            return ImageFont.truetype(path, size)
        except OSError:
            continue
    return ImageFont.load_default()

@st.cache_data
def get_mol_drawing_bytes(smiles):
    """Berechnet das RDKit-Molekül und gibt die Rohdaten zurück (lässt sich besser cachen als PIL-Bilder)."""
    if not smiles or str(smiles).lower() == "nan": 
        return None
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
    d = rdMolDraw2D.MolDraw2DCairo(500, 500) 
    opts = d.drawOptions()
    opts.bondLineWidth = 3
    opts.minFontSize = 25
    d.DrawMolecule(mol)
    d.FinishDrawing()
    return d.GetDrawingText()

def generate_mol_img(smiles):
    """Konvertiert die gecachten Rohdaten in ein PIL-Bild."""
    drawing_bytes = get_mol_drawing_bytes(smiles)
    if drawing_bytes:
        return Image.open(io.BytesIO(drawing_bytes))
    return None

def export_to_rxn(smiles_str):
    try:
        rxn = AllChem.ReactionFromSmarts(smiles_str, useSmiles=True)
        return AllChem.ReactionToRxnBlock(rxn)
    except Exception: # Generelles Error-Handling statt nacktem 'except:'
        return None

def parse_smiles(smiles_str):
    if not smiles_str or str(smiles_str).lower() == "nan": return None
    parts = str(smiles_str).split(">")
    parts = [p for p in parts if p.strip() != ""]
    
    if len(parts) == 3:
        return {"r": parts[0].split("."), "a": parts[1], "p": parts[2].split(".")}
    elif len(parts) == 2:
        return {"r": parts[0].split("."), "a": "", "p": parts[1].split(".")}
    return None

def draw_text_centered(draw, txt, x, y, font):
    """Hilfsfunktion ausgelagert für mehr Übersichtlichkeit."""
    bbox = draw.textbbox((0, 0), txt, font=font)
    th = bbox[3] - bbox[1]
    draw.text((x, y - (th // 2) - 10), txt, fill="black", font=font)
    return bbox[2] - bbox[0]

def draw_arrow(draw, x_start, y_mid, arrow_len, reagent_text=None):
    """Zeichnet den Reaktionspfeil und den Reagenz-Text."""
    draw.line([(x_start, y_mid), (x_start + arrow_len, y_mid)], fill="black", width=4)
    draw.polygon([(x_start + arrow_len, y_mid), (x_start + arrow_len - 30, y_mid - 15), (x_start + arrow_len - 30, y_mid + 15)], fill="black")
    
    if reagent_text and reagent_text.strip():
        f_reag = get_font(45)
        bbox_r = draw.textbbox((0, 0), reagent_text, font=f_reag)
        r_w = bbox_r[2] - bbox_r[0]
        draw.text((x_start + (arrow_len // 2) - (r_w // 2), y_mid - 70), reagent_text, fill="black", font=f_reag)

def draw_single_line(data, ratio_str, f_main, f_title, name=None):
    r_parts = str(ratio_str).replace(">>", ".").split(".")
    anzahl = len(data['r']) + len(data['p'])
    canvas_w = (anzahl * 600) + 1200 
    canvas_h = 850
    canvas = Image.new('RGB', (canvas_w, canvas_h), color=(255, 255, 255))
    draw = ImageDraw.Draw(canvas)
    y_mid = canvas_h // 2
    x_cursor = 100

    if name:
        draw.text((100, 30), str(name), fill="black", font=f_title)

    # --- 1. REAKTANDEN ---
    for i, sm in enumerate(data['r']):
        coeff = r_parts[i].strip() if i < len(r_parts) else ""
        if coeff not in ["1", "nan", "", "1.0"]:
            x_cursor += draw_text_centered(draw, coeff, x_cursor, y_mid, f_main) + 30
        
        m_img = generate_mol_img(sm)
        if m_img:
            canvas.paste(m_img, (x_cursor, y_mid - 250), m_img.convert("RGBA"))
            x_cursor += 500
        else:
            st.warning(f"SMILES '{sm}' in Reaktand {i+1} ist ungültig und wurde übersprungen.")
            x_cursor += 150 # Platzhalter
            
        if i < len(data['r']) - 1:
            x_cursor += 30
            x_cursor += draw_text_centered(draw, "+", x_cursor, y_mid, f_main) + 50

    # --- 2. PFEIL ---
    x_cursor += 60
    arrow_len = 350
    draw_arrow(draw, x_cursor, y_mid, arrow_len, data.get("a"))
    x_cursor += arrow_len + 70

    # --- 3. PRODUKTE ---
    offset = len(data['r'])
    for i, sm in enumerate(data['p']):
        c_idx = offset + i
        coeff = r_parts[c_idx].strip() if c_idx < len(r_parts) else ""
        if coeff not in ["1", "nan", "", "1.0"]:
            x_cursor += draw_text_centered(draw, coeff, x_cursor, y_mid, f_main) + 30
            
        m_img = generate_mol_img(sm)
        if m_img:
            canvas.paste(m_img, (x_cursor, y_mid - 250), m_img.convert("RGBA"))
            x_cursor += 500
        else:
            st.warning(f"SMILES '{sm}' in Produkt {i+1} ist ungültig und wurde übersprungen.")
            x_cursor += 150
            
        if i < len(data['p']) - 1:
            x_cursor += 30
            x_cursor += draw_text_centered(draw, "+", x_cursor, y_mid, f_main) + 50

    return canvas.crop((0, 0, x_cursor + 150, canvas_h))

def combine_and_render(smiles_list, ratio_list, reaction_name):
    f_title = get_font(60)
    f_main = get_font(90)
    images = []
    for i, sm in enumerate(smiles_list):
        data = parse_smiles(sm)
        if data:
            images.append(draw_single_line(data, ratio_list[i], f_main, f_title, reaction_name if i == 0 else None))
    if not images: return None
    max_w = max(img.width for img in images)
    total_h = sum(img.height for img in images) + 20
    combined = Image.new('RGB', (max_w, total_h), (255, 255, 255))
    y_off = 0
    for img in images:
        combined.paste(img, (0, y_off))
        y_off += img.height
    return combined

# --- 3. INTERFACE ---
st.title("⚗️ Chemie-Designer Pro")

# Neuen Tab hinzugefügt
t1, t2, t3 = st.tabs(["✨ Einzel-Eingabe", "📂 Batch-Export", "🖼️ Bild-Erkennung"])

with t1:
    name = st.text_input("Reaktionsname", "Reaktion_1")
    st.info("Format-Tipp: Edukt1.Edukt2 > Reagenz > Produkt")
    c1, c2 = st.columns(2)
    s1 = c1.text_input("SMILES 1", "C*(=O)O.OCC > H2SO4 > C*(=O)OCC")
    r1 = c1.text_input("Verhältnis 1", "1>>1")
    s2 = c2.text_input("SMILES 2 (Optional)", "")
    r2 = c2.text_input("Verhältnis 2", "1>>1")

    if st.button("Generieren"):
        s_list = [s for s in [s1, s2] if s.strip()]
        img = combine_and_render(s_list, [r1, r2], name)
        if img:
            st.image(img)
            b_png = io.BytesIO(); img.save(b_png, format="PNG")
            cd1, cd2 = st.columns(2)
            cd1.download_button("💾 Bild (PNG)", b_png.getvalue(), f"{name}.png")
            rxn = export_to_rxn(s1)
            if rxn: cd2.download_button("🧬 Editierbar (RXN)", rxn, f"{name}.rxn")

with t2:
    up = st.file_uploader("Excel hochladen", type=["xlsx"])
    if up:
        df = pd.read_excel(up)
        cols = df.columns.tolist()
        cn = st.selectbox("Name Spalte", cols); cs1 = st.selectbox("SMILES A Spalte", cols)
        cs2 = st.selectbox("SMILES B Spalte (Opt)", ["Keine"] + cols); cr = st.selectbox("Verhältnis Spalte", ["Keine"] + cols)

        if st.button("🚀 Reaktionen generieren"):
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
                for idx, row in df.iterrows():
                    nom_raw = str(row[cn])
                    nom = nom_raw.replace("/", "_").replace("\\", "_")
                    s_list = [str(row[cs1])]
                    ratio = str(row[cr]) if cr != "Keine" else "1>>1"
                    if cs2 != "Keine" and pd.notna(row[cs2]): s_list.append(str(row[cs2]))
                    
                    img = combine_and_render(s_list, [ratio, ratio], nom_raw)
                    if img:
                        st.markdown(f"### {idx+1}. {nom_raw}")
                        st.image(img)
                        img_byte_arr = io.BytesIO(); img.save(img_byte_arr, format='PNG')
                        rxn_data = export_to_rxn(str(row[cs1]))
                        c_d1, c_d2 = st.columns(2)
                        c_d1.download_button(f"PNG: {nom}", img_byte_arr.getvalue(), f"{nom}.png", key=f"p_{idx}")
                        if rxn_data:
                            c_d2.download_button(f"RXN: {nom}", rxn_data, f"{nom}.rxn", key=f"r_{idx}")
                        zip_file.writestr(f"{nom}.png", img_byte_arr.getvalue())
                        if rxn_data: zip_file.writestr(f"{nom}.rxn", rxn_data)
                        st.divider()
            
            st.success("✅ Alle Reaktionen fertig erstellt!")
            st.download_button(
                label="📦 ALLE REAKTIONEN ALS ZIP HERUNTERLADEN",
                data=zip_buffer.getvalue(),
                file_name="chemie_batch_export.zip",
                mime="application/zip",
                use_container_width=True
            )

with t3:
    st.header("Struktur aus Bild erkennen (OCSR)")
    st.info("💡 **Hinweis:** Dies erfordert die Bibliothek `DECIMER`. Wenn du sie noch nicht hast, installiere sie über das Terminal mit `pip install decimer`.")
    
    uploaded_image = st.file_uploader("Bild einer chemischen Struktur hochladen", type=["png", "jpg", "jpeg"])
    
    if uploaded_image:
        st.image(uploaded_image, caption="Hochgeladenes Bild", width=400)
        
        if st.button("🔍 Struktur erkennen"):
            try:
                # Dynamischer Import, damit die App nicht abstürzt, falls DECIMER fehlt
                from DECIMER import predict_SMILES
                
                with st.spinner("Analysiere Bild durch KI... Dies kann einen Moment dauern."):
                    # DECIMER benötigt einen Dateipfad, daher speichern wir das Bild kurz temporär
                    with tempfile.NamedTemporaryFile(delete=False, suffix=".png") as tmp:
                        tmp.write(uploaded_image.getbuffer())
                        tmp_path = tmp.name
                    
                    try:
                        predicted_smiles = predict_SMILES(tmp_path)
                        st.success("✅ Struktur erfolgreich erkannt!")
                        st.code(predicted_smiles, language="text")
                        
                        st.markdown("### Generierte Vorschau aus erkanntem SMILES")
                        redrawn_img = generate_mol_img(predicted_smiles)
                        if redrawn_img:
                            st.image(redrawn_img)
                        else:
                            st.warning("Das erkannte SMILES konnte von RDKit leider nicht gezeichnet werden.")
                    finally:
                        # Temporäre Datei wieder löschen, um Müll zu vermeiden
                        os.remove(tmp_path) 
            
            except ImportError:
                st.error("🚨 Die Bibliothek `DECIMER` wurde nicht gefunden. Bitte stoppe Streamlit, führe `pip install decimer` aus und starte die App neu.")
            except Exception as e:
                st.error(f"❌ Ein Fehler ist bei der Vorhersage aufgetreten: {e}")
