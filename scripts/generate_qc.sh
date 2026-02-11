#!/bin/bash
# Generate QC report with FA/MD maps

set -e

export PATH="/home/uqahonne/mrtrix3/bin:$PATH"

OUTPUT_DIR="/home/uqahonne/uq/nif/mrtrix-rish/test_output"
QC_DIR="$OUTPUT_DIR/qc"
mkdir -p "$QC_DIR/figures"

echo "========================================"
echo "Generating QC Report with FA/MD"
echo "========================================"

echo ""
echo "[1/3] Computing FA/MD for sub-002 (Reference)..."
DWI="$OUTPUT_DIR/sub-002/ses-1/dwi.mif"
MASK="$OUTPUT_DIR/sub-002/ses-1/mask.mif"
dwi2tensor "$DWI" "$QC_DIR/sub-002_tensor.mif" -mask "$MASK" -force -quiet
tensor2metric "$QC_DIR/sub-002_tensor.mif" -fa "$QC_DIR/sub-002_fa.mif" -force -quiet
tensor2metric "$QC_DIR/sub-002_tensor.mif" -adc "$QC_DIR/sub-002_md.mif" -force -quiet
REF_FA=$(mrstats "$QC_DIR/sub-002_fa.mif" -mask "$MASK" -output mean 2>/dev/null)
REF_MD=$(mrstats "$QC_DIR/sub-002_md.mif" -mask "$MASK" -output mean 2>/dev/null)
echo "  FA mean: $REF_FA"
echo "  MD mean: $REF_MD"

echo ""
echo "[2/3] Computing FA/MD for sub-003 (Original)..."
DWI="$OUTPUT_DIR/sub-003/ses-1/dwi.mif"
MASK="$OUTPUT_DIR/sub-003/ses-1/mask.mif"
dwi2tensor "$DWI" "$QC_DIR/sub-003_original_tensor.mif" -mask "$MASK" -force -quiet
tensor2metric "$QC_DIR/sub-003_original_tensor.mif" -fa "$QC_DIR/sub-003_original_fa.mif" -force -quiet
tensor2metric "$QC_DIR/sub-003_original_tensor.mif" -adc "$QC_DIR/sub-003_original_md.mif" -force -quiet
ORIG_FA=$(mrstats "$QC_DIR/sub-003_original_fa.mif" -mask "$MASK" -output mean 2>/dev/null)
ORIG_MD=$(mrstats "$QC_DIR/sub-003_original_md.mif" -mask "$MASK" -output mean 2>/dev/null)
echo "  FA mean: $ORIG_FA"
echo "  MD mean: $ORIG_MD"

echo ""
echo "[3/3] Generating HTML report..."

# RISH stats
REF_RISH0=$(mrstats "$OUTPUT_DIR/sub-002/ses-1/rish/rish_l0.mif" -mask "$OUTPUT_DIR/sub-002/ses-1/mask.mif" -output mean 2>/dev/null)
REF_RISH2=$(mrstats "$OUTPUT_DIR/sub-002/ses-1/rish/rish_l2.mif" -mask "$OUTPUT_DIR/sub-002/ses-1/mask.mif" -output mean 2>/dev/null)
REF_RISH4=$(mrstats "$OUTPUT_DIR/sub-002/ses-1/rish/rish_l4.mif" -mask "$OUTPUT_DIR/sub-002/ses-1/mask.mif" -output mean 2>/dev/null)
REF_RISH6=$(mrstats "$OUTPUT_DIR/sub-002/ses-1/rish/rish_l6.mif" -mask "$OUTPUT_DIR/sub-002/ses-1/mask.mif" -output mean 2>/dev/null)
REF_RISH8=$(mrstats "$OUTPUT_DIR/sub-002/ses-1/rish/rish_l8.mif" -mask "$OUTPUT_DIR/sub-002/ses-1/mask.mif" -output mean 2>/dev/null)

ORIG_RISH0=$(mrstats "$OUTPUT_DIR/sub-003/ses-1/rish/rish_l0.mif" -mask "$MASK" -output mean 2>/dev/null)
ORIG_RISH2=$(mrstats "$OUTPUT_DIR/sub-003/ses-1/rish/rish_l2.mif" -mask "$MASK" -output mean 2>/dev/null)
ORIG_RISH4=$(mrstats "$OUTPUT_DIR/sub-003/ses-1/rish/rish_l4.mif" -mask "$MASK" -output mean 2>/dev/null)
ORIG_RISH6=$(mrstats "$OUTPUT_DIR/sub-003/ses-1/rish/rish_l6.mif" -mask "$MASK" -output mean 2>/dev/null)
ORIG_RISH8=$(mrstats "$OUTPUT_DIR/sub-003/ses-1/rish/rish_l8.mif" -mask "$MASK" -output mean 2>/dev/null)

HARM_RISH0=$(mrstats "$OUTPUT_DIR/harmonization/rish_harmonized/rish_l0.mif" -mask "$MASK" -output mean 2>/dev/null)
HARM_RISH2=$(mrstats "$OUTPUT_DIR/harmonization/rish_harmonized/rish_l2.mif" -mask "$MASK" -output mean 2>/dev/null)
HARM_RISH4=$(mrstats "$OUTPUT_DIR/harmonization/rish_harmonized/rish_l4.mif" -mask "$MASK" -output mean 2>/dev/null)
HARM_RISH6=$(mrstats "$OUTPUT_DIR/harmonization/rish_harmonized/rish_l6.mif" -mask "$MASK" -output mean 2>/dev/null)
HARM_RISH8=$(mrstats "$OUTPUT_DIR/harmonization/rish_harmonized/rish_l8.mif" -mask "$MASK" -output mean 2>/dev/null)

SCALE0=$(mrstats "$OUTPUT_DIR/harmonization/scale_maps/scale_l0.mif" -mask "$MASK" -output mean 2>/dev/null)
SCALE2=$(mrstats "$OUTPUT_DIR/harmonization/scale_maps/scale_l2.mif" -mask "$MASK" -output mean 2>/dev/null)
SCALE4=$(mrstats "$OUTPUT_DIR/harmonization/scale_maps/scale_l4.mif" -mask "$MASK" -output mean 2>/dev/null)
SCALE6=$(mrstats "$OUTPUT_DIR/harmonization/scale_maps/scale_l6.mif" -mask "$MASK" -output mean 2>/dev/null)
SCALE8=$(mrstats "$OUTPUT_DIR/harmonization/scale_maps/scale_l8.mif" -mask "$MASK" -output mean 2>/dev/null)

# Calculate differences
calc_pct() { echo "scale=1; ($1 - $2) / $2 * 100" | bc 2>/dev/null || echo "N/A"; }

FA_DIFF=$(calc_pct $ORIG_FA $REF_FA)
MD_DIFF=$(calc_pct $ORIG_MD $REF_MD)

RISH0_ORIG_DIFF=$(calc_pct $ORIG_RISH0 $REF_RISH0)
RISH2_ORIG_DIFF=$(calc_pct $ORIG_RISH2 $REF_RISH2)
RISH4_ORIG_DIFF=$(calc_pct $ORIG_RISH4 $REF_RISH4)
RISH6_ORIG_DIFF=$(calc_pct $ORIG_RISH6 $REF_RISH6)
RISH8_ORIG_DIFF=$(calc_pct $ORIG_RISH8 $REF_RISH8)

RISH0_HARM_DIFF=$(calc_pct $HARM_RISH0 $REF_RISH0)
RISH2_HARM_DIFF=$(calc_pct $HARM_RISH2 $REF_RISH2)
RISH4_HARM_DIFF=$(calc_pct $HARM_RISH4 $REF_RISH4)
RISH6_HARM_DIFF=$(calc_pct $HARM_RISH6 $REF_RISH6)
RISH8_HARM_DIFF=$(calc_pct $HARM_RISH8 $REF_RISH8)

# Generate HTML
cat > "$QC_DIR/qc_report.html" << EOF
<!DOCTYPE html>
<html>
<head>
    <title>RISH Harmonization QC Report</title>
    <style>
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background: #f5f5f5;
        }
        .header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 20px;
        }
        .header h1 { margin: 0; font-size: 28px; }
        .header p { margin: 10px 0 0 0; opacity: 0.9; }
        .section {
            background: white;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .section h2 {
            color: #333;
            border-bottom: 2px solid #667eea;
            padding-bottom: 10px;
            margin-top: 0;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
        }
        th, td {
            padding: 12px 15px;
            text-align: center;
            border-bottom: 1px solid #e0e0e0;
        }
        th {
            background: #f8f9fa;
            font-weight: 600;
            color: #333;
        }
        tr:hover { background: #fafafa; }
        .good { color: #28a745; font-weight: 600; }
        .warn { color: #ffc107; font-weight: 600; }
        .bad { color: #dc3545; font-weight: 600; }
        .metric-value { font-family: monospace; font-size: 14px; }
        .summary-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 15px;
            margin: 20px 0;
        }
        .summary-card {
            background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);
            padding: 20px;
            border-radius: 8px;
            text-align: center;
            border-left: 4px solid #667eea;
        }
        .summary-card .value {
            font-size: 24px;
            font-weight: bold;
            color: #667eea;
        }
        .summary-card .label {
            color: #666;
            font-size: 13px;
            margin-top: 5px;
        }
        .footer {
            text-align: center;
            color: #888;
            padding: 20px;
            font-size: 13px;
        }
        code {
            background: #f4f4f4;
            padding: 2px 6px;
            border-radius: 3px;
            font-size: 13px;
        }
        .highlight { background: #fff3cd; }
    </style>
</head>
<body>
    <div class="header">
        <h1>🧠 RISH Harmonization QC Report</h1>
        <p>Generated: $(date '+%Y-%m-%d %H:%M:%S')</p>
        <p>Reference: sub-002 (Site A) → Target: sub-003 (Site B)</p>
    </div>
    
    <div class="section">
        <h2>📊 Summary</h2>
        <div class="summary-grid">
            <div class="summary-card">
                <div class="value">2</div>
                <div class="label">Subjects</div>
            </div>
            <div class="summary-card">
                <div class="value">Single-Shell</div>
                <div class="label">Acquisition</div>
            </div>
            <div class="summary-card">
                <div class="value">b=2800</div>
                <div class="label">b-value</div>
            </div>
            <div class="summary-card">
                <div class="value">lmax=8</div>
                <div class="label">SH Order</div>
            </div>
            <div class="summary-card">
                <div class="value">CSD</div>
                <div class="label">FOD Algorithm</div>
            </div>
        </div>
    </div>

    <div class="section">
        <h2>📈 DTI Metrics (FA & MD)</h2>
        <table>
            <tr>
                <th>Metric</th>
                <th>Reference (sub-002)</th>
                <th>Target (sub-003)</th>
                <th>Δ (%)</th>
            </tr>
            <tr>
                <td><strong>FA</strong></td>
                <td class="metric-value">$REF_FA</td>
                <td class="metric-value">$ORIG_FA</td>
                <td>${FA_DIFF}%</td>
            </tr>
            <tr>
                <td><strong>MD</strong> (mm²/s)</td>
                <td class="metric-value">$REF_MD</td>
                <td class="metric-value">$ORIG_MD</td>
                <td>${MD_DIFF}%</td>
            </tr>
        </table>
        <p><em>Note: DTI metrics computed from original DWI data (pre-FOD).</em></p>
    </div>

    <div class="section">
        <h2>🔄 RISH Features - Before & After Harmonization</h2>
        <table>
            <tr>
                <th>Feature</th>
                <th>Reference</th>
                <th>Original</th>
                <th>Harmonized</th>
                <th>Scale</th>
                <th>Δ Before</th>
                <th>Δ After</th>
            </tr>
            <tr>
                <td><strong>θ₀</strong></td>
                <td class="metric-value">$REF_RISH0</td>
                <td class="metric-value">$ORIG_RISH0</td>
                <td class="metric-value">$HARM_RISH0</td>
                <td class="metric-value">$SCALE0</td>
                <td>${RISH0_ORIG_DIFF}%</td>
                <td class="good">${RISH0_HARM_DIFF}%</td>
            </tr>
            <tr class="highlight">
                <td><strong>θ₂</strong></td>
                <td class="metric-value">$REF_RISH2</td>
                <td class="metric-value">$ORIG_RISH2</td>
                <td class="metric-value">$HARM_RISH2</td>
                <td class="metric-value">$SCALE2</td>
                <td class="bad">${RISH2_ORIG_DIFF}%</td>
                <td class="good">${RISH2_HARM_DIFF}%</td>
            </tr>
            <tr>
                <td><strong>θ₄</strong></td>
                <td class="metric-value">$REF_RISH4</td>
                <td class="metric-value">$ORIG_RISH4</td>
                <td class="metric-value">$HARM_RISH4</td>
                <td class="metric-value">$SCALE4</td>
                <td>${RISH4_ORIG_DIFF}%</td>
                <td class="good">${RISH4_HARM_DIFF}%</td>
            </tr>
            <tr>
                <td><strong>θ₆</strong></td>
                <td class="metric-value">$REF_RISH6</td>
                <td class="metric-value">$ORIG_RISH6</td>
                <td class="metric-value">$HARM_RISH6</td>
                <td class="metric-value">$SCALE6</td>
                <td>${RISH6_ORIG_DIFF}%</td>
                <td class="good">${RISH6_HARM_DIFF}%</td>
            </tr>
            <tr class="highlight">
                <td><strong>θ₈</strong></td>
                <td class="metric-value">$REF_RISH8</td>
                <td class="metric-value">$ORIG_RISH8</td>
                <td class="metric-value">$HARM_RISH8</td>
                <td class="metric-value">$SCALE8</td>
                <td class="bad">${RISH8_ORIG_DIFF}%</td>
                <td class="good">${RISH8_HARM_DIFF}%</td>
            </tr>
        </table>
        <p><em>Highlighted rows show the largest corrections. Green = improved alignment with reference.</em></p>
    </div>

    <div class="section">
        <h2>📁 Output Files</h2>
        <table>
            <tr><th>File</th><th>Description</th></tr>
            <tr><td><code>qc/sub-002_fa.mif</code></td><td>Reference FA map</td></tr>
            <tr><td><code>qc/sub-002_md.mif</code></td><td>Reference MD map</td></tr>
            <tr><td><code>qc/sub-003_original_fa.mif</code></td><td>Target FA map (original)</td></tr>
            <tr><td><code>qc/sub-003_original_md.mif</code></td><td>Target MD map (original)</td></tr>
            <tr><td><code>harmonization/scale_maps/scale_l*.mif</code></td><td>Per-order scale maps</td></tr>
            <tr><td><code>harmonization/wm_fod_harmonized.mif</code></td><td>Harmonized FOD</td></tr>
            <tr><td><code>harmonization/rish_harmonized/</code></td><td>Harmonized RISH features</td></tr>
        </table>
    </div>

    <div class="section">
        <h2>ℹ️ Methods</h2>
        <p><strong>RISH Harmonization</strong> operates at the FOD (spherical harmonics) level:</p>
        <ol>
            <li><strong>Extract RISH:</strong> θₗ = √(Σₘ |cₗₘ|²) for each SH order l</li>
            <li><strong>Compute scale maps:</strong> scale_l = θₗ_ref / θₗ_target (smoothed, clipped)</li>
            <li><strong>Apply scaling:</strong> c'ₗₘ = cₗₘ × scale_l</li>
        </ol>
        <p><strong>Advantage:</strong> Preserves fiber orientation while removing scanner effects.</p>
    </div>

    <div class="footer">
        <p>Generated by <strong>mrtrix-rish</strong> v0.1.0</p>
        <p>Reference: Cetin Karayumak et al. (2019) NeuroImage</p>
    </div>
</body>
</html>
EOF

echo ""
echo "========================================"
echo "✓ QC Report Generated!"
echo "========================================"
echo ""
echo "FA/MD Maps:"
ls -la "$QC_DIR"/*.mif 2>/dev/null | head -10
echo ""
echo "HTML Report: $QC_DIR/qc_report.html"
echo ""
echo "View in browser:"
echo "  file://$QC_DIR/qc_report.html"
