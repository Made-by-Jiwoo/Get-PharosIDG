import streamlit as st
import pandas as pd
import requests

# --- Pharos API 호출 함수 ---
def get_pharos_data(gene_symbols):
    url = 'https://pharos-api.ncats.io/graphql'
    results = {}
    query = """
    query getTarget($gene: String!) {
      target(q: { sym: $gene }) {
        sym, name, tdl, fam, uniprot
        publications(top: 10) { pmid, title, journal, date }
      }
    }
    """
    for gene in gene_symbols:
        variables = {'gene': gene}
        try:
            response = requests.post(url, json={'query': query, 'variables': variables}, timeout=10)
            if response.status_code == 200:
                resp_json = response.json()
                target_data = resp_json.get('data', {}).get('target')
                results[gene] = target_data if target_data else {'error': 'No data found'}
            else:
                results[gene] = {'error': f'HTTP {response.status_code}'}
        except Exception as e:
            results[gene] = {'error': str(e)}
    return results

# --- Open Targets API 호출 함수 ---
def get_opentargets_data(uniprot_id):
    if not uniprot_id: return None
    url = "https://api.platform.opentargets.org/api/v4/graphql"
    query = """
    query targetByUniprot($uId: [String!]!) {
      mapIds(queryTerms: $uId) {
        mappings { hits { object { ... on Target { id, approvedSymbol
                knownDrugs { count, rows { drug { name }, phase, status } }
                associatedDiseases { count } } } } } }
    }
    """
    variables = {"uId": [uniprot_id]}
    try:
        response = requests.post(url, json={'query': query, 'variables': variables}, timeout=10)
        if response.status_code == 200:
            res_json = response.json()
            mappings = res_json.get('data', {}).get('mapIds', {}).get('mappings', [])
            if mappings and mappings[0].get('hits'):
                return mappings[0]['hits'][0].get('object')
    except: pass
    return None

# --- AlphaFold API 호출 함수 ---
def get_alphafold_pdb(uniprot_id):
    if not uniprot_id: return None, None
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        response = requests.get(api_url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            if data and len(data) > 0:
                pdb_url = data[0].get('pdbUrl')
                pdb_content = requests.get(pdb_url).content
                return pdb_url, pdb_content
    except: pass
    return None, None

# --- Streamlit UI 설정 ---
st.set_page_config(page_title="Biobytes Target Analyzer", layout="wide")

st.title("Pharos외에 AlphaFold PDB를 가져오기 ")
st.markdown("유전자명을 입력하면 **IDG Level**을 확인하고 **AlphaFold의 DB에서 PDB**를 즉시 다운로드 ㄱㄱ.")

# 세션 상태 초기화
if 'analysis_results' not in st.session_state:
    st.session_state.analysis_results = None

input_text = st.text_input("유전자 기호 입력 (예: ETS2, EGFR) 쉼표로 구분을 하니까 여러개도 한버넹 쓰세요. 추가로 띄어쓰기까지는 포함이 가능", placeholder="ETS2, EGFR")

# 버튼 클릭 시에만 API 호출 후 세션 상태에 저장
if st.button("데이터 분석 및 PDB 찾기"):
    if input_text:
        gene_list = [g.strip().upper() for g in input_text.split(",") if g.strip()]
        with st.spinner('데이터 분석 중...'):
            # API 데이터를 가져와서 세션에 저장
            st.session_state.analysis_results = get_pharos_data(gene_list)
    else:
        st.warning("유전자 기호를 입력해 주세요.")

# 세션에 데이터가 있을 경우 화면에 표시 (버튼 클릭 여부와 상관없이 유지됨)
if st.session_state.analysis_results:
    pharos_info = st.session_state.analysis_results
    summary_data = []
    valid_genes = []

    for gene, info in pharos_info.items():
        if info and 'error' not in info:
            summary_data.append({
                "Gene Symbol": gene,
                "Full Name": info.get('name'),
                "TDL": info.get('tdl'),
                "Family": info.get('fam'),
                "UniProt ID": info.get('uniprot')
            })
            valid_genes.append(gene)
        else:
            st.error(f"**{gene}**: 데이터를 찾을 수 없습니다.")

    if summary_data:
        st.subheader("📊 분석 요약 결과")
        st.dataframe(pd.DataFrame(summary_data), use_container_width=True)
        st.divider()

        # 이제 여기서 다른 유전자를 선택해도 데이터가 사라지지 않습니다.
        selected_gene = st.selectbox("상세 리포트를 볼 유전자를 선택하세요", valid_genes)

        if selected_gene:
            info = pharos_info[selected_gene]
            uniprot_id = info.get('uniprot')
            
            # 상세 정보 로딩 (캐싱 가능)
            ot_data = get_opentargets_data(uniprot_id)
            pdb_url, pdb_content = get_alphafold_pdb(uniprot_id)

            st.markdown(f"### 🧬 {selected_gene} 통합 리포트")
            
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Full Name", info.get('name', 'N/A'))
                st.metric("Family", info.get('fam', 'N/A'))
            with col2:
                st.metric("TDL (개발단계)", info.get('tdl', 'N/A'), help="Tbio: 생물학적 연구 위주, Tchem: 화합물 존재")
                st.write(f"**UniProt ID**: `{uniprot_id}`")
            with col3:
                if ot_data:
                    st.metric("질병 연관성", f"{ot_data.get('associatedDiseases', {}).get('count', 0)} 건")
                    st.metric("알려진 약물", f"{ot_data.get('knownDrugs', {}).get('count', 0)} 건")

            tab1, tab2, tab3 = st.tabs(["📚 최근 관련 논문", "💊 약물 현황", "🔬 AlphaFold PDB"])
            
            with tab1:
                st.subheader("최근 관련 논문 (Top 10)")
                pubs = info.get('publications', [])
                if pubs:
                    for pub in pubs:
                        date_str = str(pub['date'])[:4] if pub.get('date') else "N/A"
                        st.markdown(f"- **({date_str})** {pub['title']}  \n  *Journal: {pub['journal']}* | [PMID: {pub['pmid']}](https://pubmed.ncbi.nlm.nih.gov/{pub['pmid']}/)")
                else: st.info("관련 논문 정보가 없습니다.")

            with tab2:
                st.subheader("약물 및 임상 현황")
                if ot_data and ot_data.get('knownDrugs', {}).get('count', 0) > 0:
                    drugs = ot_data['knownDrugs']['rows']
                    drug_df = pd.DataFrame([
                        {"Drug Name": d['drug']['name'], "Phase": d['phase'], "Status": d['status']}
                        for d in drugs[:10]
                    ])
                    st.table(drug_df)
                else: st.info("알려진 약물 정보가 없습니다.")

            with tab3:
                st.subheader("AlphaFold 구조 데이터")
                if pdb_url:
                    st.success(f"PDB 파일을 찾았습니다: [링크]({pdb_url})")
                    st.download_button(
                        label=f"{selected_gene} PDB 다운로드",
                        data=pdb_content,
                        file_name=f"AF_{selected_gene}_{uniprot_id}.pdb",
                        mime="application/octet-stream",
                        key=f"dl_{selected_gene}"
                    )
                else: st.error("AlphaFold PDB 정보를 찾을 수 없습니다.")

st.divider()
st.caption("Integrated by Biobytes | Data from Pharos, Open Targets & AlphaFold DB")
