import streamlit as st
import requests
import pandas as pd
import io

def get_pharos_data(gene_symbols):
    url = 'https://pharos-api.ncats.io/graphql'
    results = {}

    # 1. 쿼리 수정: q 대신 직접 sym으로 검색하거나 facet 활용
    # 2. properties 대신 더 직접적인 명칭 사용
    query = """
    query getTarget($gene: String!) {
      target(q: { sym: $gene }) {
        sym
        tdl
        uniprot  # Pharos 스키마에서 uniprot은 보통 직접 호출 가능합니다.
      }
    }
    """

    for gene in gene_symbols:
        variables = {'gene': gene}
        try:
            response = requests.post(url, json={'query': query, 'variables': variables}, timeout=10)
            
            if response.status_code == 200:
                resp_json = response.json()
                
                # 에러 메시지가 있는지 확인
                if 'errors' in resp_json:
                    print(f"Error for {gene}: {resp_json['errors'][0]['message']}")
                    continue
                
                target_data = resp_json.get('data', {}).get('target')
                if target_data:
                    results[gene] = {
                        'tdl': target_data.get('tdl'),
                        'uniprot': target_data.get('uniprot')
                    }
                else:
                    print(f"No data found for {gene}")
            else:
                print(f"HTTP Error {response.status_code} for {gene}")
        except Exception as e:
            print(f"Exception for {gene}: {str(e)}")
            continue
            
    return results

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

st.title("🧬 AlphaFold-Pharos Integrated Analyzer")
st.markdown("유전자명을 입력하면 **IDG Level**을 확인하고 **AlphaFold PDB**를 즉시 다운로드합니다.")

input_text = st.text_input("유전자 기호 입력 (예: ETS2, EGFR)", placeholder="ETS2, EGFR")

if st.button("데이터 분석 및 PDB 찾기"):
    if input_text:
        gene_list = [g.strip().upper() for g in input_text.split(",") if g.strip()]
        
        with st.spinner('데이터를 통합 분석 중입니다...'):
            pharos_info = get_pharos_data(gene_list)
            
            final_results = []
            
            for gene in gene_list:
                info = pharos_info.get(gene, {'tdl': 'Not Found', 'uniprot': None})
                tdl = info['tdl']
                uniprot_id = info['uniprot']
                
                pdb_url, pdb_content = get_alphafold_pdb(uniprot_id)
                
                final_results.append({
                    "Gene": gene,
                    "IDG Level": tdl,
                    "UniProt ID": uniprot_id if uniprot_id else "N/A",
                    "AF PDB Link": pdb_url if pdb_url else "Not Found",
                    "pdb_content": pdb_content
                })

            # 테이블 표시
            df = pd.DataFrame(final_results).drop(columns=['pdb_content'])
            st.subheader("📊 분석 결과 요약")
            st.table(df)

            # 다운로드 섹션
            st.subheader("📥 PDB 구조 파일 다운로드")
            cols = st.columns(len(final_results))
            
            for idx, item in enumerate(final_results):
                with cols[idx]:
                    st.write(f"**{item['Gene']}**")
                    if item['pdb_content']:
                        st.download_button(
                            label=f"Download PDB",
                            data=item['pdb_content'],
                            file_name=f"AF_{item['Gene']}_{item['UniProt ID']}.pdb",
                            mime="application/octet-stream",
                            key=f"btn_{idx}"
                        )
                    else:
                        st.error("PDB 없음")
    else:
        st.warning("유전자 기호를 입력해 주세요.")

st.divider()
st.caption("Integrated by Biobytes | Data from Pharos & AlphaFold DB (EBI)")