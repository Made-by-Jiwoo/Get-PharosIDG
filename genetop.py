import streamlit as st
import pandas as pd
import requests

# --- Pharos API 호출 함수 ---
def get_pharos_target_levels(gene_symbols):
    url = 'https://pharos-api.ncats.io/graphql'
    results = []

    query = """
    query getTarget($gene: String!) {
      target(q: { sym: $gene }) {
        sym
        tdl
      }
    }
    """

    for gene in gene_symbols:
        variables = {'gene': gene}
        try:
            response = requests.post(url, json={'query': query, 'variables': variables}, timeout=10)
            if response.status_code == 200:
                data = response.json()
                target_data = data.get('data', {}).get('target')
                if target_data:
                    results.append(target_data)
                else:
                    results.append({'sym': gene, 'tdl': 'Not Found'})
            else:
                results.append({'sym': gene, 'tdl': f'Error ({response.status_code})'})
        except Exception as e:
            results.append({'sym': gene, 'tdl': f'Connection Error'})

    return results

def get_pharos_data(gene_symbols):
    url = 'https://pharos-api.ncats.io/graphql'
    results = {}

    query = """
    query getTarget($gene: String!) {
      target(q: { sym: $gene }) {
        sym
        name
        tdl
        fam
        uniprot
        publications(top: 10) {
          pmid
          title
          journal
          date
        }
      }
    }
    """

    for gene in gene_symbols:
        variables = {'gene': gene}
        try:
            response = requests.post(url, json={'query': query, 'variables': variables}, timeout=10)
            if response.status_code == 200:
                resp_json = response.json()
                if 'errors' in resp_json:
                    results[gene] = {'error': resp_json['errors'][0]['message']}
                    continue
                
                target_data = resp_json.get('data', {}).get('target')
                if target_data:
                    results[gene] = target_data
                else:
                    results[gene] = {'error': 'No data found'}
            else:
                results[gene] = {'error': f'HTTP {response.status_code}'}
        except Exception as e:
            results[gene] = {'error': str(e)}
            
    return results

# --- Open Targets API 호출 함수 ---
def get_opentargets_data(uniprot_id):
    if not uniprot_id:
        return None
    
    url = "https://api.platform.opentargets.org/api/v4/graphql"
    query = """
    query targetByUniprot($uId: [String!]!) {
      mapIds(queryTerms: $uId) {
        mappings {
          hits {
            object {
              ... on Target {
                id
                approvedSymbol
                knownDrugs {
                  count
                  rows {
                    drug { name }
                    phase
                    status
                  }
                }
                associatedDiseases {
                  count
                }
              }
            }
          }
        }
      }
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
    except:
        pass
    return None

# --- 2. AlphaFold API: PDB 파일 링크 및 데이터 가져오기 ---
def get_alphafold_pdb(uniprot_id):
    if not uniprot_id:
        return None, None
    
    # AlphaFold API 호출
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        response = requests.get(api_url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            if data and len(data) > 0:
                pdb_url = data[0].get('pdbUrl')
                # 실제 PDB 파일 내용 다운로드
                pdb_content = requests.get(pdb_url).content
                return pdb_url, pdb_content
    except:
        pass
    return None, None

# --- Streamlit UI 설정 ---
st.set_page_config(page_title="Biobytes Target Analyzer", layout="wide")

st.title("Pharos외에 AlphaFold PDB를 가져오기 ")
st.markdown("유전자명을 입력하면 **IDG Level**을 확인하고 **AlphaFold의 DB에서 PDB**를 즉시 다운로드 ㄱㄱ.")

input_text = st.text_input("유전자 기호 입력 (예: ETS2, EGFR) 쉼표로 구분을 하니까 여러개도 한버넹 쓰세요. 추가로 띄어쓰기까지는 포함이 가능", placeholder="ETS2, EGFR")

if st.button("데이터 분석 및 PDB 찾기"):
    if input_text:
        gene_list = [g.strip().upper() for g in input_text.split(",") if g.strip()]
        
        with st.spinner('데이터 분석 중...'):
            pharos_info = get_pharos_data(gene_list)
            
            for gene in gene_list:
                info = pharos_info.get(gene, {})
                
                if 'error' in info:
                    st.error(f"**{gene}**: {info['error']}")
                    continue

                uniprot_id = info.get('uniprot')
                ot_data = get_opentargets_data(uniprot_id)
                pdb_url, pdb_content = get_alphafold_pdb(uniprot_id)

                # 토글(Expander)로 감싸기
                with st.expander(f"🧬 {gene} 통합 리포트 (클릭하여 상세 보기)", expanded=False):
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

                    # 상세 정보 탭
                    tab1, tab2, tab3 = st.tabs(["📚 최근 관련 논문", "💊 약물 현황", "🔬 AlphaFold PDB"])
                    
                    with tab1:
                        st.subheader("최근 관련 논문 (Top 10)")
                        pubs = info.get('publications', [])
                        if pubs:
                            for pub in pubs:
                                date_str = str(pub['date'])[:4] if pub.get('date') else "N/A"
                                st.markdown(f"- **({date_str})** {pub['title']}  \n  *Journal: {pub['journal']}* | [PMID: {pub['pmid']}](https://pubmed.ncbi.nlm.nih.gov/{pub['pmid']}/)")
                        else:
                            st.info("관련 논문 정보가 없습니다.")

                    with tab2:
                        st.subheader("약물 및 임상 현황")
                        if ot_data and ot_data.get('knownDrugs', {}).get('count', 0) > 0:
                            drugs = ot_data['knownDrugs']['rows']
                            drug_df = pd.DataFrame([
                                {"Drug Name": d['drug']['name'], "Phase": d['phase'], "Status": d['status']}
                                for d in drugs[:10]
                            ])
                            st.table(drug_df)
                        else:
                            st.info("알려진 약물 정보가 없습니다.")

                    with tab3:
                        st.subheader("AlphaFold 구조 데이터")
                        if pdb_url:
                            st.success(f"PDB 파일을 찾았습니다: [링크]({pdb_url})")
                            st.download_button(
                                label=f"{gene} PDB 다운로드",
                                data=pdb_content,
                                file_name=f"AF_{gene}_{uniprot_id}.pdb",
                                mime="application/octet-stream",
                                key=f"dl_{gene}"
                            )
                        else:
                            st.error("AlphaFold PDB 정보를 찾을 수 없습니다.")

    else:
        st.warning("유전자 기호를 입력해 주세요.")

st.divider()
st.caption("Integrated by Biobytes | Data from Pharos, Open Targets & AlphaFold DB")
